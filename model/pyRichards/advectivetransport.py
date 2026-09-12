import numpy as np
import numba
from collections import Counter
from mpi4py import MPI

@numba.njit #Avoid parallel=True because the advection loop updates shared donor/receiver masses.
def _advect_tracer_step_kernel(tracer_mass, q_m3s, water_volume, dt, track_pair_mass):
    nrisfu, nsoil = tracer_mass.shape
    mass = tracer_mass.copy()

    if track_pair_mass:
        pair_mass = np.zeros((nrisfu, nrisfu, nsoil), dtype=np.float64)
    else:
        pair_mass = np.empty((0, 0, 0), dtype=np.float64)

    for layer in range(nsoil):
        for donor_candidate in range(nrisfu):
            for receiver_candidate in range(donor_candidate + 1, nrisfu):
                flow = q_m3s[donor_candidate, receiver_candidate, layer]
                if abs(flow) < 1e-20:
                    continue

                if flow > 0.0:
                    donor = donor_candidate
                    receiver = receiver_candidate
                else:
                    donor = receiver_candidate
                    receiver = donor_candidate

                requested_volume = abs(flow) * dt
                available_volume = max(water_volume[donor, layer], 0.0)
                volume = min(requested_volume, available_volume)
                if volume <= 0.0:
                    continue

                donor_concentration = mass[donor, layer] / water_volume[donor, layer]
                moved_mass = min(volume * donor_concentration, mass[donor, layer])

                mass[donor, layer] -= moved_mass
                mass[receiver, layer] += moved_mass
                water_volume[donor, layer] -= volume
                water_volume[receiver, layer] += volume

                if track_pair_mass:
                    pair_mass[donor, receiver, layer] += moved_mass

    return mass, pair_mass


@numba.njit
def _redistribute_mass_hrus_kernel(mass_hrus, delta_risfu, farea):
    nrisfu, nhrus = farea.shape
    nsoil = mass_hrus.shape[1]
    mass_hrus_updated = mass_hrus.copy()

    for risfu_index in range(nrisfu):
        for tracer_index in range(nsoil):
            delta = delta_risfu[risfu_index, tracer_index]
            if abs(delta) <= 1e-12:
                continue

            total_weight = 0.0
            if delta > 0.0:
                for hru_index in range(nhrus):
                    if farea[risfu_index, hru_index] > 0.0:
                        total_weight += farea[risfu_index, hru_index]
            else:
                for hru_index in range(nhrus):
                    if (farea[risfu_index, hru_index] > 0.0 and
                            abs(mass_hrus[hru_index, tracer_index]) > 1e-12):
                        total_weight += mass_hrus[hru_index, tracer_index]

                if total_weight <= 1e-12:
                    raise ValueError("RISFU has a mass loss but no HRU with mass")

            for hru_index in range(nhrus):
                if delta > 0.0:
                    if farea[risfu_index, hru_index] > 0.0:
                        weight = farea[risfu_index, hru_index] / total_weight
                        mass_hrus_updated[hru_index, tracer_index] += delta * weight
                elif (farea[risfu_index, hru_index] > 0.0 and
                      abs(mass_hrus[hru_index, tracer_index]) > 1e-12):
                    weight = mass_hrus[hru_index, tracer_index] / total_weight
                    mass_hrus_updated[hru_index, tracer_index] += delta * weight

    return mass_hrus_updated

class AdvectiveTransport:
    """
    Conservative tracer advection using richards q links.
    """
    def __init__(self, nhrus, nrisfus, nsoil, reg_ids, comm):
        self.nhrus = nhrus
        self.nrisfus = nrisfus
        self.nsoil = nsoil
        self.c0_hrus = np.zeros((nhrus, nsoil), dtype=float)
        self.tracer_mass_risfu = np.zeros((nrisfus, nsoil), dtype=float)
        self.tracer_mass_int_new = np.zeros((nrisfus, nsoil), dtype=float)
        self.tracer_mass_reg_new = np.zeros((nrisfus, nsoil), dtype=float)
        self.tracer_mass_hrus_new = np.zeros((nhrus, nsoil), dtype=float)
        self.tracer_mass_hrus_old = np.zeros((nhrus, nsoil), dtype=float)
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
        self.tracer_mass_reg_risfu = np.zeros((n,self.nsoil))

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

    def aggregate_mass_risfu(self, tracer_mass_hrus, farea):
        """
        Aggregate HRU tracer mass to RISFU tracer mass.

        `farea` is the HRU-to-RISFU area-fraction matrix with shape
        (nrisfu, nhrus).
        """
        tracer_mass_hrus = np.asarray(tracer_mass_hrus, dtype=float)
        farea = np.asarray(farea, dtype=float)

        return (farea > 0.0).astype(int) @ tracer_mass_hrus

    def compute_int_tracer(self, tracer_mass_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt, store_mass_transfer=False):
        """
        Compute the intermediate tracer mass in each RISFU.
        
        When `store_mass_transfer=True`, save a donor->receiver matrix for this step.
        """
        if store_mass_transfer:
            tracer_mass_risfu_new, pair_mass = self.advect_tracer_step(tracer_mass_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt, track_pair_mass=True)
            self.int_mass_transfer_matrix = pair_mass
        else:
            tracer_mass_risfu_new = self.advect_tracer_step(tracer_mass_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt)
            self.int_mass_transfer_matrix = None

        return tracer_mass_risfu_new

    def compute_reg_tracer(self, tracer_mass_risfu, q_m3s, local_risfu_indices, area_risfu, theta_risfu, dz_risfu, dt, store_mass_transfer=False):
        """
        Compute the regional tracer mass in each RISFU.

        `local_risfu_indices` identifies the subset of RISFUs that belong to the
        local region; return only those rows in the original regional ordering.
        
        When `store_mass_transfer=True`, save the full global pair_mass matrix.
        """
        if store_mass_transfer:
            tracer_mass_risfu_reg_new, pair_mass = self.advect_tracer_step(tracer_mass_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt, track_pair_mass=True)
            self.reg_mass_transfer_matrix = pair_mass
        else:
            tracer_mass_risfu_reg_new = self.advect_tracer_step(tracer_mass_risfu, q_m3s, area_risfu, theta_risfu, dz_risfu, dt)
            self.reg_mass_transfer_matrix = None
            
        local_risfu_indices = np.asarray(local_risfu_indices, dtype=int)
        local_tracer_mass = np.zeros((len(local_risfu_indices), self.nsoil), dtype=float)
        for i, idx in enumerate(local_risfu_indices):
            local_tracer_mass[i, :] = tracer_mass_risfu_reg_new[idx, :]

        return local_tracer_mass

    def compute_loc_tracer(self, tracer_mass_hrus_hbands, q_links, area_hrus_hbands, theta_hrus_hbands, dz_hrus_hbands, dt, store_mass_transfer=False):
        """
        Compute local tracer advection.

        When `store_mass_transfer=True`, save a donor->receiver matrix for this
        step where entry [i, j] is tracer mass moved from HRU i to HRU j.
        """
        if store_mass_transfer:
            tracer_mass_hrus_hbands_new, pair_mass = self.advect_tracer_step(
                tracer_mass_hrus_hbands,
                q_links,
                area_hrus_hbands,
                theta_hrus_hbands,
                dz_hrus_hbands,
                dt,
                track_pair_mass=True,
            )

            if tracer_mass_hrus_hbands.shape[0] == self.nhrus:
                self.local_mass_transfer_matrix = pair_mass
                self.local_mass_divergence[:] = np.sum(pair_mass, axis=1)
                self.local_mass_convergence[:] = np.sum(pair_mass, axis=0)
            else:
                # Keep diagnostics disabled for non-HRU local schemes (e.g., hbands).
                self.clear_local_mass_transfer()
        else:
            tracer_mass_hrus_hbands_new = self.advect_tracer_step(tracer_mass_hrus_hbands, q_links, area_hrus_hbands, theta_hrus_hbands, dz_hrus_hbands, dt)

        return tracer_mass_hrus_hbands_new

    def redistribute_mass_hrus(self, mass_hrus, new_mass_risfu, old_mass_risfu, farea):
        """
        Redistribute RISFU tracer mass back to HRUs.

        `farea` is the HRU-to-RISFU area-fraction matrix with shape
        (nrisfu, nhrus).
        """
        new_mass_risfu = np.asarray(new_mass_risfu, dtype=float)
        old_mass_risfu = np.asarray(old_mass_risfu, dtype=float)
        farea = np.asarray(farea, dtype=float)
        mass_hrus = np.asarray(mass_hrus, dtype=float)

        delta_risfu = new_mass_risfu - old_mass_risfu
        mass_hrus_updated = _redistribute_mass_hrus_kernel(
            mass_hrus,
            delta_risfu,
            farea,
        )

        # Check that the updated HRU masses still aggregate to the requested RISFU masses.
        assert np.all(mass_hrus_updated >= -1e-10)
        assert np.allclose((farea > 0.0).astype(int) @ mass_hrus_updated, new_mass_risfu)
        return mass_hrus_updated #farea.T @ tracer_mass_risfu

    def advect_tracer_step(self, tracer_mass, q_m3s, area_risfu, theta_risfu, dz_risfu, dt, track_pair_mass=False):
        """
        Conservatively advect tracer mass for one step using Richards q links.

        `tracer_mass` is mass per node and soil layer. The returned array uses
        the same units and shape.
        """
        tracer_mass = np.asarray(tracer_mass, dtype=float)
        area_risfu = np.asarray(area_risfu, dtype=float)
        theta_risfu = np.asarray(theta_risfu, dtype=float)
        dz_risfu = np.asarray(dz_risfu, dtype=float)
        q_m3s = np.asarray(q_m3s, dtype=float)
        nrisfu, nsoil = tracer_mass.shape

        if theta_risfu.shape != (nrisfu, nsoil):
            raise ValueError("theta must have shape (nrisfu, nsoil).")
        if dz_risfu.shape != (nrisfu, nsoil):
            raise ValueError("dz must have shape (nrisfu, nsoil).")
        if area_risfu.shape != (nrisfu,):
            raise ValueError("area must have shape (nrisfu,).")

        # Water volume per node/layer [m3]
        water_volume = area_risfu[:, None] * dz_risfu * theta_risfu
        water_volume = np.maximum(water_volume, 1e-20)

        conc_new, pair_mass = _advect_tracer_step_kernel(
            tracer_mass,
            q_m3s,
            water_volume,
            dt,
            track_pair_mass,
        )
        if track_pair_mass:
            return conc_new, pair_mass
        return conc_new

def exchange_tracer_mass_regional_units(cids, rank, HBdb):
    # This function exchanges tracer mass between regional units across ranks.
    # It uses non-blocking MPI communication to send and receive tracer mass based on the connections defined in risfu_mapping.
    if len(cids) == 0:
        return

    cid_rank_mapping = HBdb[cids[0]].cid_rank_mapping
    cids_core = set(cids)
    comm = HBdb[cids[0]].advectivetransport.comm

    # Keep tracer-mass exchanges in an exclusive tag space to avoid cross-variable collisions.
    mass_tag_offset = 700000

    request_send = []
    request_recv = []
    mass_db_all = {}

    # Track communication plans so mismatches can fail fast instead of hanging in Waitall.
    planned_send = []
    planned_recv = []

    for cid_in_rank in cids:
        local_subsurface = HBdb[cid_in_rank].advectivetransport
        tracer_mass_risfu = local_subsurface.tracer_mass_risfu[:]
        risfu_mapping = HBdb[cid_in_rank].mssubsurface.risfu_mapping

        mass_db = {cid_in_rank: tracer_mass_risfu}

        for cid in risfu_mapping.keys():
            if cid not in cids_core:
                if cid_in_rank in risfu_mapping[cid].keys():
                    rows_to_send = risfu_mapping[cid][cid_in_rank]
                    data_to_send = tracer_mass_risfu[rows_to_send, :]
                    dest = cid_rank_mapping[cid]
                    tag = mass_tag_offset + (cid_in_rank * 1000 + cid)

                    request_send.append(local_subsurface.comm.Isend(data_to_send, dest=dest, tag=tag))
                    planned_send.append((rank, dest, tag, data_to_send.shape[0], data_to_send.shape[1], str(data_to_send.dtype)))

            elif cid != cid_in_rank and cid in risfu_mapping[cid_in_rank].keys():
                rows = risfu_mapping[cid_in_rank][cid]
                mass_other_cid = HBdb[cid].advectivetransport.tracer_mass_risfu[:]
                mass_db[cid] = mass_other_cid[rows, :]

        for cid in risfu_mapping.keys():
            if cid not in cids_core and cid in risfu_mapping[cid_in_rank].keys():
                src = cid_rank_mapping[cid]
                tag = mass_tag_offset + (cid * 1000 + cid_in_rank)
                rows = risfu_mapping[cid_in_rank][cid]
                # Match dtype with local mass arrays to avoid MPI datatype inconsistencies.
                data_recv_buffer = np.empty((len(rows), tracer_mass_risfu.shape[1]), dtype=tracer_mass_risfu.dtype)

                request_recv.append(local_subsurface.comm.Irecv(data_recv_buffer, source=src, tag=tag))
                planned_recv.append((src, rank, tag, data_recv_buffer.shape[0], data_recv_buffer.shape[1], str(data_recv_buffer.dtype)))
                mass_db[cid] = data_recv_buffer

            mass_db_all[cid_in_rank] = mass_db

    gathered_sends = comm.allgather(planned_send)
    gathered_recvs = comm.allgather(planned_recv)
    global_send_plan = sum(gathered_sends, [])
    global_recv_plan = sum(gathered_recvs, [])
    if Counter(global_send_plan) != Counter(global_recv_plan):
        missing_recvs = list((Counter(global_send_plan) - Counter(global_recv_plan)).elements())[:5]
        missing_sends = list((Counter(global_recv_plan) - Counter(global_send_plan)).elements())[:5]
        raise ValueError(
            f"Tracer-mass MPI plan mismatch on rank={rank}. "
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
        mass_db = mass_db_all[cid_in_rank]

        tracer_mass = np.empty((1, local_subsurface.nsoil))
        for cid in local_subsurface.reg_ids.keys():
            if cid not in mass_db:
                raise ValueError(f"Missing tracer mass data for cid_in_rank={cid_in_rank}, cid={cid} on rank={rank}")
            tracer_mass = np.concatenate((tracer_mass, mass_db[cid]), axis=0)

        local_subsurface.tracer_mass_reg_risfu[:, :] = tracer_mass[1:, :]

    return