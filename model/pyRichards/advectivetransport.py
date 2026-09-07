import numpy as np
import numba
from collections import Counter
from mpi4py import MPI

class AdvectiveTransport:
    """
    Conservative tracer advection using richards q links.
    """
    def __init__(self, nhrus, nrisfus, nsoil, reg_ids, comm):
        self.nhrus = nhrus
        self.nrisfus = nrisfus
        self.nsoil = nsoil
        self.c0_hrus = np.zeros((nhrus, nsoil), dtype=float)
        self.c_risfu = np.zeros((nrisfus, nsoil), dtype=float)
        self.c_int_new = np.zeros((nrisfus, nsoil), dtype=float)
        self.c_reg_new = np.zeros((nrisfus, nsoil), dtype=float)
        self.c_hrus_new = np.zeros((nhrus, nsoil), dtype=float)
        self.c_hbands = None
        # Per-step local HRU donor->receiver tracer mass transfer [mass/time step].
        # Filled only when the HRU local-flow scheme requests tracking.
        self.local_mass_transfer_matrix = None
        self.local_mass_convergence = np.zeros((nhrus, nsoil), dtype=float)
        self.local_mass_divergence = np.zeros((nhrus, nsoil), dtype=float)
        # Per-step intermediate RISFU donor->receiver mass transfer [mass/time step].
        self.int_mass_transfer_matrix = None
        # Per-step regional RISFU donor->receiver mass transfer [mass/time step].
        self.reg_mass_transfer_matrix = None
        self.comm = comm
        self.reg_ids = reg_ids

        #Initialize regional matrices
        n = 0
        for cid in self.reg_ids.keys():
            n += len(self.reg_ids[cid])
        self.reg_conc_risfu = np.zeros((n,self.nsoil))

    def clear_local_mass_transfer(self):
        """
        Reset per-step local HRU mass-transfer diagnostics.
        """
        self.local_mass_transfer_matrix = None
        self.local_mass_convergence[:] = 0.0
        self.local_mass_divergence[:] = 0.0

    def synthetic_concentration(self,):
        """
        Build a synthetic initial concentration field [mass / m3 water].
        Higher near highest HRU index, decays with depth.
        """
        x = np.arange(self.nhrus, dtype=float)
        center = 0.9 * max(1, self.nhrus - 1)
        sigma = max(1.0, 0.25 * self.nhrus)

        horizontal = np.exp(-0.5 * ((x - center) / sigma) ** 2)
        horizontal /= max(horizontal.max(), 1e-20)

        depth_scale = np.linspace(1.0, 0.7, self.nsoil)
        c0 = horizontal[:, None] * depth_scale[None, :]

        return c0

    def aggregate_concentration_risfu(self, c_hrus, farea, area_hrus=None, theta_hrus=None, dz_hrus=None):
        """
        Aggregate HRU concentrations to RISFU concentrations using the mass-conservative
        water volume in each HRU layer. `farea` is the area-fraction matrix with shape
        (nrisfu, nhrus).

        If the water-volume inputs are not supplied, fall back to the legacy area-only
        weighted average for compatibility with older call sites.
        """
        c_hrus = np.asarray(c_hrus, dtype=float)
        farea = np.asarray(farea, dtype=float)

        if area_hrus is None or theta_hrus is None or dz_hrus is None:
            return farea @ c_hrus

        area_hrus = np.asarray(area_hrus, dtype=float)
        theta_hrus = np.asarray(theta_hrus, dtype=float)
        dz_hrus = np.asarray(dz_hrus, dtype=float)

        water_hrus = area_hrus[:, None] * theta_hrus * dz_hrus
        mass_hrus = c_hrus * water_hrus

        mass_risfu = farea @ mass_hrus
        water_risfu = farea @ water_hrus

        return mass_risfu / np.maximum(water_risfu, 1e-20)

    def compute_int_tracer(self, c_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt, store_mass_transfer=False):
        """
        Compute the intermediate tracer mass in each RISFU.
        
        When `store_mass_transfer=True`, save a donor->receiver matrix for this step.
        """
        if store_mass_transfer:
            c_risfu_new, pair_mass = self.advect_tracer_step(c_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt, track_pair_mass=True)
            self.int_mass_transfer_matrix = pair_mass
        else:
            c_risfu_new = self.advect_tracer_step(c_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt)
            self.int_mass_transfer_matrix = None

        return c_risfu_new

    def compute_reg_tracer(self, c_risfu, q_m3s, local_risfu_indices, area_risfu, theta_risfu, dz_risfu, dt, store_mass_transfer=False):
        """
        Compute the regional tracer mass in each RISFU.

        `local_risfu_indices` identifies the subset of RISFUs that belong to the
        local region; return only those rows in the original regional ordering.
        
        When `store_mass_transfer=True`, save the full global pair_mass matrix.
        """
        if store_mass_transfer:
            c_risfu_reg_new, pair_mass = self.advect_tracer_step(c_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt, track_pair_mass=True)
            self.reg_mass_transfer_matrix = pair_mass
        else:
            c_risfu_reg_new = self.advect_tracer_step(c_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt)
            self.reg_mass_transfer_matrix = None
            
        local_risfu_indices = np.asarray(local_risfu_indices, dtype=int)
        local_net_tracer_conc = np.zeros((len(local_risfu_indices), self.nsoil), dtype=float)
        for i, idx in enumerate(local_risfu_indices):
            local_net_tracer_conc[i, :] = c_risfu_reg_new[idx, :]

        return local_net_tracer_conc

    def compute_loc_tracer(self, c_hrus_hbands, q_links, area_hrus_hbands, theta_hrus_hbands, dz_hrus_hbands, dt, store_mass_transfer=False):
        """
        Compute local tracer advection.

        When `store_mass_transfer=True`, save a donor->receiver matrix for this
        step where entry [i, j] is tracer mass moved from HRU i to HRU j.
        """
        if store_mass_transfer:
            c_hbands_hrus_new, pair_mass = self.advect_tracer_step(
                c_hrus_hbands,
                q_links,
                area_hrus_hbands,
                theta_hrus_hbands,
                dz_hrus_hbands,
                dt,
                track_pair_mass=True,
            )

            if c_hrus_hbands.shape[0] == self.nhrus:
                self.local_mass_transfer_matrix = pair_mass
                self.local_mass_divergence[:] = np.sum(pair_mass, axis=1)
                self.local_mass_convergence[:] = np.sum(pair_mass, axis=0)
            else:
                # Keep diagnostics disabled for non-HRU local schemes (e.g., hbands).
                self.clear_local_mass_transfer()
        else:
            c_hbands_hrus_new = self.advect_tracer_step(c_hrus_hbands, q_links, area_hrus_hbands, theta_hrus_hbands, dz_hrus_hbands, dt)

        return c_hbands_hrus_new

    def redistribute_concentration_hrus(self, c_risfu, farea, area_risfu=None, theta_risfu=None, dz_risfu=None, area_hrus=None, theta_hrus=None, dz_hrus=None):
        """
        Redistribute RISFU concentrations back to HRUs using the water volume in each
        RISFU layer. `farea` is the area-fraction matrix with shape (nrisfu, nhrus).

        If the water-volume inputs are not supplied, fall back to the legacy area-only
        weighting for compatibility with older call sites.
        """
        c_risfu = np.asarray(c_risfu, dtype=float)
        farea = np.asarray(farea, dtype=float)

        if area_risfu is None or theta_risfu is None or dz_risfu is None:
            return farea.T @ c_risfu

        area_risfu = np.asarray(area_risfu, dtype=float)
        theta_risfu = np.asarray(theta_risfu, dtype=float)
        dz_risfu = np.asarray(dz_risfu, dtype=float)

        water_risfu = area_risfu[:, None] * theta_risfu * dz_risfu
        mass_risfu = c_risfu * water_risfu

        mass_hrus = farea.T @ mass_risfu

        if area_hrus is None or theta_hrus is None or dz_hrus is None:
            return mass_hrus / np.maximum(np.sum(water_risfu, axis=0, keepdims=True).T, 1e-20)

        area_hrus = np.asarray(area_hrus, dtype=float)
        theta_hrus = np.asarray(theta_hrus, dtype=float)
        dz_hrus = np.asarray(dz_hrus, dtype=float)

        water_hrus = area_hrus[:, None] * theta_hrus * dz_hrus
        return mass_hrus / np.maximum(water_hrus, 1e-20)

    def advect_tracer_step(self, c_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt, track_pair_mass=False):
        """
        Conservative tracer advection for one step using richards q links.
        """
        c_risfu = np.asarray(c_risfu, dtype=float)
        area_risfu = np.asarray(area_risfu, dtype=float)
        theta_risfu = np.asarray(theta_risfu, dtype=float)
        dz_risfu = np.asarray(dz_risfu, dtype=float)
        q_m3s = np.asarray(q_m3s, dtype=float)
        nrisfu, nsoil = c_risfu.shape

        if theta_risfu.shape != (nrisfu, nsoil):
            raise ValueError("theta must have shape (nrisfu, nsoil).")
        if dz_risfu.shape != (nrisfu, nsoil):
            raise ValueError("dz must have shape (nrisfu, nsoil).")
        if area_risfu.shape != (nrisfu,):
            raise ValueError("area must have shape (nrisfu,).")

        # Water volume per node/layer [m3]
        Vw = area_risfu[:, None] * dz_risfu * theta_risfu
        Vw = np.maximum(Vw, 1e-20)

        # Tracer mass [mass]
        M = c_risfu #* Vw #Previous mass
        pair_mass = None
        if track_pair_mass:
            pair_mass = np.zeros((nrisfu, nrisfu, nsoil), dtype=float)

        for il in range(nsoil):
            for i in range(nrisfu):
                for j in range(i + 1, nrisfu):
                    f = q_m3s[i, j, il]  # [m3/s] positive means i -> j
                    if abs(f) < 1e-20:
                        continue

                    if f > 0.0:
                        donor = i
                        recv = j
                    else:
                        donor = j
                        recv = i

                    qabs = abs(f)
                    vol = qabs * dt

                    c_donor = M[donor, il] / Vw[donor, il]
                    m_potential = vol * c_donor
                    m_move = min(m_potential, M[donor, il])

                    M[donor, il] -= m_move
                    M[recv, il] += m_move
                    Vw[donor, il] -= vol
                    Vw[recv, il] += vol

                    if track_pair_mass:
                        pair_mass[donor, recv, il] += m_move

        conc_new = M #/ Vw
        if track_pair_mass:
            return conc_new, pair_mass
        return conc_new

def exchange_concentrations_regional_units(cids, rank, HBdb):
    # This function performs the exchange of tracer concentrations between regional units across ranks.
    # It uses non-blocking MPI communication to send and receive the necessary concentrations based on the connections defined in the risfu_mapping.
    if len(cids) == 0:
        return

    cid_rank_mapping = HBdb[cids[0]].cid_rank_mapping
    cids_core = set(cids)
    comm = HBdb[cids[0]].advectivetransport.comm

    # Keep concentrations exchanges in an exclusive tag space to avoid cross-variable collisions.
    conc_tag_offset = 700000

    request_send = []
    request_recv = []
    conc_db_all = {}

    # Track communication plans so mismatches can fail fast instead of hanging in Waitall.
    planned_send = []
    planned_recv = []

    for cid_in_rank in cids:
        local_subsurface = HBdb[cid_in_rank].advectivetransport
        conc_hb = local_subsurface.c_risfu[:]
        risfu_mapping = HBdb[cid_in_rank].mssubsurface.risfu_mapping

        conc_db = {cid_in_rank: conc_hb}

        for cid in risfu_mapping.keys():
            if cid not in cids_core:
                if cid_in_rank in risfu_mapping[cid].keys():
                    rows_to_send = risfu_mapping[cid][cid_in_rank]
                    data_to_send = conc_hb[rows_to_send, :]
                    dest = cid_rank_mapping[cid]
                    tag = conc_tag_offset + (cid_in_rank * 1000 + cid)

                    request_send.append(local_subsurface.comm.Isend(data_to_send, dest=dest, tag=tag))
                    planned_send.append((rank, dest, tag, data_to_send.shape[0], data_to_send.shape[1], str(data_to_send.dtype)))

            elif cid != cid_in_rank and cid in risfu_mapping[cid_in_rank].keys():
                rows = risfu_mapping[cid_in_rank][cid]
                conc_other_cid = HBdb[cid].advectivetransport.c_risfu[:]
                conc_db[cid] = conc_other_cid[rows, :]

        for cid in risfu_mapping.keys():
            if cid not in cids_core and cid in risfu_mapping[cid_in_rank].keys():
                src = cid_rank_mapping[cid]
                tag = conc_tag_offset + (cid * 1000 + cid_in_rank)
                rows = risfu_mapping[cid_in_rank][cid]
                # Match dtype with local conc arrays to avoid MPI datatype inconsistencies.
                data_recv_buffer = np.empty((len(rows), conc_hb.shape[1]), dtype=conc_hb.dtype)

                request_recv.append(local_subsurface.comm.Irecv(data_recv_buffer, source=src, tag=tag))
                planned_recv.append((src, rank, tag, data_recv_buffer.shape[0], data_recv_buffer.shape[1], str(data_recv_buffer.dtype)))
                conc_db[cid] = data_recv_buffer

        conc_db_all[cid_in_rank] = conc_db

    gathered_sends = comm.allgather(planned_send)
    gathered_recvs = comm.allgather(planned_recv)
    global_send_plan = sum(gathered_sends, [])
    global_recv_plan = sum(gathered_recvs, [])
    if Counter(global_send_plan) != Counter(global_recv_plan):
        missing_recvs = list((Counter(global_send_plan) - Counter(global_recv_plan)).elements())[:5]
        missing_sends = list((Counter(global_recv_plan) - Counter(global_send_plan)).elements())[:5]
        raise ValueError(
            f"Concentration MPI plan mismatch on rank={rank}. "
            f"Missing recv matches (sample): {missing_recvs}; "
            f"missing send matches (sample): {missing_sends}"
        )

    comm.Barrier()
    all_requests = request_recv + request_send
    if len(all_requests) > 0:
        MPI.Request.Waitall(all_requests)
    comm.Barrier()

    for cid_in_rank in cids:
        local_subsurface = HBdb[cid_in_rank].advectivetransport
        conc_db = conc_db_all[cid_in_rank]

        conc = np.empty((1, local_subsurface.nsoil))
        for cid in local_subsurface.reg_ids.keys():
            if cid not in conc_db:
                raise ValueError(f"Missing concentration data for cid_in_rank={cid_in_rank}, cid={cid} on rank={rank}")
            conc = np.concatenate((conc, conc_db[cid]), axis=0)

        local_subsurface.reg_conc_risfu[:, :] = conc[1:, :]

    return