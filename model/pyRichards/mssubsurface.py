import warnings
warnings.filterwarnings('ignore')
from mpi4py import MPI
import numpy as np
import numba
import pickle
from collections import Counter

class mssubsurface:
    def __init__(self,info,comm,current_cid):
        #Get the communicator
        self.cid  = current_cid #call the cid   
        self.comm = comm
        #Open the file with the area fraction for intercells
        farea           = pickle.load(open(info['gw_area_file'],'rb'))
        self.total_area = np.sum(farea, axis = 1)
        self.units      = np.count_nonzero(np.sum(farea, axis = 1))
        self.farea_gw   = np.zeros((self.units,farea.shape[1]))
        self.farea_gw[:]= farea[self.total_area > 0]
        #Open files with characteristics of groundwater units
        self.area_units = pickle.load(open(info['gw_areg_file'],'rb')) 
        self.w_gw       = pickle.load(open(info['gw_wreg_file'],'rb'))
        self.dx_gw      = pickle.load(open(info['gw_dxrg_file'],'rb'))       
        self.mconx      = pickle.load(open(info['gw_conx_file'],'rb')) #mconx intermediate level
        self.dem_gw     = pickle.load(open(info['gw_elev_file'],'rb'))        
        #Open files with characteristics of groundwater units for macroscale polygons
        self.ccid       = pickle.load(open(info['gw_ccid_file'],'rb')) #mconx regional level
        self.dxcid      = pickle.load(open(info['gw_dxcd_file'],'rb')) #dx_cids
        self.wcid       = pickle.load(open(info['gw_wcid_file'],'rb')) #w_cids
        self.reg_ids    = pickle.load(open(info['gw_rids_file'],'rb'))
        #Open the parameters file and average them for the regional/intermediate unit. Parameters are for current cid
        dz              = np.array(info['dz'])
        self.nsoil      = dz.size
        self.dz_gw      = np.zeros((self.units,self.nsoil)) + dz
        param           = pickle.load(open(info['gw_parm_file'],'rb'))
        self.tr_gw      = param['tr'][:,:] #residual soil moisture
        self.ts_gw      = param['ts'][:,:] #saturated soil moisture
        self.bb_gw      = param['bb'][:,:] #curve parameter Brooks and Corey
        self.ks_gw      = param['ks'][:,:] #saturated hydraulic conductivity
        self.sp_gw      = param['sp'][:,:] #saturated matric potential
        self.m_gw       = param['m'] #curve parameter Brooks and Corey
        
        # Extract saturation flux flag from JSON configuration
        self.flag_sat = info.get('multiscale_subsurface', {}).get('flag_saturatedflux', False)
       
        #Define parameters for regional interaction
        self.subdomain_interaction_parameters(info)
        
        #Initialize regional matrices
        n = 0
        for cid in self.reg_ids.keys():
            n += len(self.reg_ids[cid])
        self.reg_theta_gw = np.zeros((n,self.nsoil))
        self.reg_temperature_gw = np.zeros((n,self.nsoil))
        self.reg_dz_gw = np.zeros((n,self.nsoil)) + dz

        # Initialize multiscale components
        self.hdiv_loc = None
        self.hdiv_int = None
        self.hdiv_reg = None
        self.hdiv_reg_same_cid = None
        self.hdiv_total = None
        self.hdiv_heat_loc = None
        self.hdiv_heat_int = None
        self.hdiv_heat_reg = None
        self.hdiv_heat_total = None
        self.q_int = None
        self.q_reg = None
        self.q_reg_same_cid = None
        self.inter_unit_flow_m3s = None
        self.regional_inter_unit_flow_m3s = None
        self.regional_inter_unit_flow_m3s_cross = None
        self.regional_inter_unit_flow_m3s_same = None
        self.reg_owner_cid = None
        self.this_cid = None

        # Regional indexing maps each exchanged unit to its source CID.
        self.reg_owner_cid, self.this_cid = self._build_regional_indexing()
        
        return
    
    def aggregate_variable(self,variable):
        aggregated_variable = parameters_aggregation(self.nsoil, variable, self.farea_gw, self.units);
        return aggregated_variable

    def _unit_divergence_from_flow_tensor(self, flow_tensor):
        # flow_tensor[from_unit, to_unit, layer] stores signed link flow in m3/s.
        # Summing over outgoing links returns unit divergence [m3/s].
        return np.sum(flow_tensor, axis=1)

    def _distribute_unit_divergence_to_hrus(self, unit_divergence_m3s, area_hrus):
        # Map unit-level divergence [m3/s] to HRU-level divergence [mm/s].
        q_hru_m3s = np.dot(self.farea_gw.T, unit_divergence_m3s)
        hdiv_hru = np.zeros_like(q_hru_m3s)
        area = area_hrus[:, np.newaxis]
        np.divide(-1000.0 * q_hru_m3s, area, out=hdiv_hru, where=area > 0.0)
        return hdiv_hru

    def _build_regional_indexing(self):
        # Build stable row-index -> CID ownership for exchanged regional units.
        reg_owner_cid = np.zeros(self.nreg_un, dtype=np.int64)
        this_cid = None
        n = 0
        for cid in self.reg_ids.keys():
            n_units = self.reg_ids[cid].size
            reg_owner_cid[n:n + n_units] = cid
            if cid == self.cid:
                this_cid = np.arange(n_units, dtype=np.int64) + n
            n += n_units
        if this_cid is None:
            raise ValueError(f"Current cid={self.cid} not found in reg_ids for regional indexing")
        return reg_owner_cid, this_cid
   
    def subdomain_interaction_parameters(self,info):
        # Load the parameters for the regional/intermediate units interacting with the current subdomain. 
        # The parameters are averaged for the regional/intermediate unit. 
        # The parameters are stored in self.reg_tr_gw, self.reg_ts_gw, self.reg_bb_gw, self.reg_ks_gw, self.reg_sp_gw, self.reg_m_gw, self.reg_dem_gw.

        #Define experiment directory
        rdir = '%s/experiments/simulations/%s' % (info['rdir'],info['experiment'])
        #Define number of regional/intermediate units and ids
        self.nreg_un = self.ccid.shape[0] #number of regional units for regional interaction
        cids         = self.reg_ids.keys() #this contains a vector with the cids interacting with the current macroscale polygon
        #Define the parameters for the set of groundwater units
        self.reg_m_gw = []; #for parameter m regional
        for cid in cids:
            #open the parameters file
            parameters = pickle.load(open('%s/%s/groundwater/param_aggregation.pck' % (rdir, cid),'rb'))
            #units of that cid interacting with the subdomain
            units_cid  = self.reg_ids[cid].astype(int) - 1; #-1 to match python indexing
            #extract the parameters for those units
            param_cid  = parameters['m'][units_cid]
            self.reg_m_gw  = np.append(self.reg_m_gw,param_cid)
            
         #for parameter dem
        self.reg_dem_gw = [];
        for cid in cids:
            #open the dem file
            dem         = pickle.load(open('%s/%s/groundwater/elv.pck' % (rdir,cid),'rb'))
            #units of that cid interacting with the subdomain
            units_cid   = self.reg_ids[cid].astype(int) - 1; #-1 to match python indexing
            #extract the parameters for those units
            param_cid   = dem[units_cid]
            self.reg_dem_gw = np.append(self.reg_dem_gw,param_cid)
            
         #for the rest of the hydraulic properties
        vars = ['tr','ts','bb','ks','sp']
        for var in vars:
            gw_param = np.empty((1,self.nsoil));
            exec('self.reg_%s_gw = np.empty((self.nreg_un,self.nsoil))' % var);
            for cid in cids:
                #open the parameters file
                parameters = pickle.load(open('%s/%s/groundwater/param_aggregation.pck' % (rdir,cid),'rb'))
                #units of that cid interacting with the subdomain
                units_cid  = self.reg_ids[cid].astype(int) - 1; #-1 to match python indexing
                #extract the parameters for those units
                param_cid  = parameters[var][units_cid,:]
                gw_param   = np.concatenate((gw_param, param_cid), axis = 0)
            exec('self.reg_%s_gw[:,:] = gw_param[1:,:]' % var);
        
        return
    
    #-------------------Compute Flows between RISFU Units---------------------------
    # Regional interaction: compute the divergence between regional units
    def compute_regional_hdiv(self, area_hrus=None):
        
        #Call the function to compute the divergence
        regional_inter_unit_flow_m3s = np.zeros((self.nreg_un, self.nreg_un, self.nsoil))
        self.regional_inter_unit_flow_m3s = update_workhorse_int_gw(self.reg_theta_gw,self.reg_dz_gw,regional_inter_unit_flow_m3s,\
                                                                    self.reg_tr_gw,self.reg_ts_gw,self.reg_bb_gw,\
                                                                    self.reg_sp_gw,self.reg_m_gw,self.reg_ks_gw,\
                                                                    self.reg_dem_gw,self.wcid,self.dxcid,\
                                                                        self.ccid,self.af,self.flag_sat)
        
        # Refresh indexing if dimensions changed.
        reg_owner_cid = getattr(self, 'reg_owner_cid', None)
        if reg_owner_cid is None or reg_owner_cid.size != self.nreg_un:
            self.reg_owner_cid, self.this_cid = self._build_regional_indexing()

        # Separate cross-CID (regional) links from same-CID links.
        cross_mask = self.reg_owner_cid[:, np.newaxis] != self.reg_owner_cid[np.newaxis, :]
        self.regional_inter_unit_flow_m3s_cross = np.where(
            cross_mask[:, :, np.newaxis],
            self.regional_inter_unit_flow_m3s,
            0.0,
        )
        self.regional_inter_unit_flow_m3s_same = np.where(
            cross_mask[:, :, np.newaxis],
            0.0,
            self.regional_inter_unit_flow_m3s,
        )

        # Keep only cross-CID divergence in the regional bucket.
        q_regional_all = self._unit_divergence_from_flow_tensor(self.regional_inter_unit_flow_m3s_cross)
        q_same_all = self._unit_divergence_from_flow_tensor(self.regional_inter_unit_flow_m3s_same)
        self.q_reg = np.copy(q_regional_all[self.this_cid, :])
        self.q_reg_same_cid = np.copy(q_same_all[self.this_cid, :])

        if area_hrus is not None:
            self.hdiv_reg = self._distribute_unit_divergence_to_hrus(self.q_reg, area_hrus)
            self.hdiv_reg_same_cid = self._distribute_unit_divergence_to_hrus(self.q_reg_same_cid, area_hrus)
        
        return

    # Backward-compatible alias used by existing tests/workflows.
    def update_subdomains_regional(self):
        return self.compute_regional_hdiv()

    def compute_regional_hdiv_heat(self, area_hrus, dz_hrus):
        flows_reg = getattr(self, 'regional_inter_unit_flow_m3s_cross', None)
        self.hdiv_heat_reg = compute_enthalpy_flux_regional(flows_reg,self.reg_temperature_gw,self.this_cid,self.farea_gw,\
                                                            area_hrus,dz_hrus)
        return np.copy(self.hdiv_heat_reg)
    
    # Intermediate interaction: compute the divergence between intermediate units within the same macroscale polygon.
    def compute_intermediate_hdiv(self, area_hrus=None):
        theta_gw  = self.th_gw; 
        af        = self.af;
        # Compute unit-to-unit flow tensor first and then distribute to HRUs.
        inter_unit_flow_m3s  = np.zeros((theta_gw.shape[0], theta_gw.shape[0], self.nsoil))
        self.inter_unit_flow_m3s = update_workhorse_int_gw(theta_gw,self.dz_gw,inter_unit_flow_m3s,self.tr_gw,self.ts_gw,\
                                                           self.bb_gw,self.sp_gw,self.m_gw,self.ks_gw,self.dem_gw,self.w_gw,self.dx_gw,\
                                                               self.mconx,af,self.flag_sat)
        self.q_int = self._unit_divergence_from_flow_tensor(self.inter_unit_flow_m3s)
        self.hdiv_int = self._distribute_unit_divergence_to_hrus(self.q_int, area_hrus)
        return

    # Backward-compatible alias used by existing tests/workflows.
    def update_subdomains_intermediate(self, area_hrus):
        return self.compute_intermediate_hdiv(area_hrus)

    def compute_intermediate_hdiv_heat(self, area_hrus, dz_hrus):
        self.hdiv_heat_int = compute_enthalpy_flux(self.farea_gw, self.inter_unit_flow_m3s,
                                                   self.temp_gw, area_hrus, dz_hrus)
        return np.copy(self.hdiv_heat_int)

#-------------------Defining the functions---------------------------
#solve richards for the intermediate interaction
@numba.jit(nopython=True, cache=True)
def update_workhorse_int_gw(theta_gw, dz_gw, inter_unit_flow_m3s, thetar_gw, thetas_gw,
                            b_gw, satpsi_gw, m_gw, ksat_gw, dem_gw, w_gw, dx_gw,
                            mconx, af, flag_sat=False):
    # flag_sat passed as parameter; Dupuit-Forchheimer approximation when True
    for il in range(theta_gw.shape[1]):
        # Calculate soil moisture potential
        psi = calculate_soil_moisture_potential(il, theta_gw, thetar_gw[:, il],thetas_gw[:, il], b_gw[:, il], satpsi_gw[:, il])
        zbot = np.sum(dz_gw[:,0:il+1],axis=1)
        ztop = zbot - dz_gw[:, il]
        T = calculate_transmissivity(psi, ztop, zbot, m_gw, ksat_gw[:, il],satpsi_gw[:, il], b_gw[:, il], af)
        h = calculate_hydraulic_head(dem_gw, psi, ztop)
        q = calculate_lateralflow_gw(h, T, w_gw, dx_gw, mconx)
        #Apply sat/unsaturated separation
        if flag_sat:
            eps = 0.01
            for i in range(q.shape[0]):
                for j in range(q.shape[1]):
                    to_unit = mconx[i,j]
                    # Suppress link if IETHER endpoint is unsaturated
                    if (theta_gw[i, il] <= (1 - eps) * thetas_gw[i, il] or 
                        theta_gw[to_unit,il] <= (1-eps)*thetas_gw[to_unit,il]):
                        q[i, j] = 0.0 #no horizontal flow in unsaturated layers
                        reverse_slots = np.where(mconx[to_unit, :] == i)[0]
                        for reverse_j in reverse_slots:
                            q[to_unit, reverse_j] = 0.0 # Also zero reverse link
        for i in range(q.shape[0]):
            for j in range(q.shape[1]):
                to_unit = mconx[i, j]
                inter_unit_flow_m3s[i, to_unit, il] += q[i, j]

    return inter_unit_flow_m3s

#solve enthalpy flux for the intermediate interaction
def compute_enthalpy_flux(clusters, flows, temperatures, area_hrus, dz_hrus):
    rho_w = 1000.0
    c_w = 4186.0
    eps = 1e-20

    nclusters, nhrus = clusters.shape
    nsoil = temperatures.shape[1]

    hdiv_heat = np.zeros((nhrus, nsoil))
    for il in range(nsoil):
        # Heat transport per cluster [J/s], sign aligned with flow direction
        net_heat_cluster = np.zeros(nclusters)
        #print(f'[intermediate flux] Layer {il}: flows = {flows[:,:,il].sum()} {np.abs(flows[:,:,il]).mean():.2e}', flush=True)
        for k in range(nclusters):
            for m in range(nclusters):
                #q_link = 0.5 * (flows[k, m, il] - flows[m, k, il])  # m3/s
                q_link = flows[k, m, il]  # Use signed flow directly for clarity
                if np.abs(q_link) <= eps:
                    continue

                if q_link > 0.0:
                    upwind_cluster = k
                    downwind_cluster = m
                else:
                    upwind_cluster = m
                    downwind_cluster = k

                transported_heat = np.abs(q_link) * rho_w * c_w * temperatures[upwind_cluster, il]
                net_heat_cluster[upwind_cluster] += transported_heat
                net_heat_cluster[downwind_cluster] -= transported_heat

        #print(f'[intermediate enthalpy] Layer {il}: net_heat_cluster = {net_heat_cluster.sum()} {np.abs(net_heat_cluster).mean():.2e}', flush=True)
        # Map cluster heat power [J/s] to HRUs [J/s]
        net_heat_hru = np.zeros(nhrus)
        for i in range(nhrus):
            for k in range(nclusters):
                net_heat_hru[i] += clusters[k, i] * net_heat_cluster[k]

        # Convert to volumetric divergence [J/m3/s]
        for i in range(nhrus):
            volume = area_hrus[i] * dz_hrus[i, il]
            if volume > eps:
                #hdiv_heat[i, il] = net_heat_hru[i] / volume
                hdiv_heat[i, il] = net_heat_hru[i] / area_hrus[i]
            else:
                hdiv_heat[i, il] = 0.0

    return hdiv_heat

#solve enthalpy flux for the regional interaction
def compute_enthalpy_flux_regional(flows, unit_temperatures, local_unit_indices, local_clusters, area_hrus, dz_hrus):
    rho_w = 1000.0
    c_w = 4186.0
    eps = 1e-20

    nclusters = flows.shape[0]
    nhrus = local_clusters.shape[1]
    nsoil = flows.shape[2]

    hdiv_heat = np.zeros((nhrus, nsoil))
    for il in range(nsoil):
        #print(f'[regional flux] Layer {il}: flows = {flows[:,:,il].sum()} {np.abs(flows[:,:,il]).mean():.2e}', flush=True)
        net_heat_cluster = np.zeros(nclusters)
        for k in range(nclusters):
            for m in range(nclusters):
                #q_link = 0.5 * (flows[k, m, il] - flows[m, k, il]) #m3/s, signed flow from k to m
                q_link = flows[k, m, il]  # Use signed flow directly for clarity
                if np.abs(q_link) <= eps:
                    continue

                if q_link > 0.0:
                    upwind_cluster = k
                    downwind_cluster = m
                else:
                    upwind_cluster = m
                    downwind_cluster = k

                transported_heat = np.abs(q_link) * rho_w * c_w * unit_temperatures[upwind_cluster, il]
                net_heat_cluster[upwind_cluster] += transported_heat
                net_heat_cluster[downwind_cluster] -= transported_heat

        #print(f'[regional enthalpy] Layer {il}: net_heat_cluster = {net_heat_cluster.sum()} {np.abs(net_heat_cluster).mean():.2e}', flush=True)
        local_net_heat_units = np.zeros(local_clusters.shape[0])
        for i in range(local_clusters.shape[0]):
            local_net_heat_units[i] = net_heat_cluster[local_unit_indices[i]]

        net_heat_hru = np.zeros(nhrus)
        for hru in range(nhrus):
            for i in range(local_clusters.shape[0]):
                net_heat_hru[hru] += local_clusters[i, hru] * local_net_heat_units[i]

        for hru in range(nhrus):
            volume = area_hrus[hru] * dz_hrus[hru, il]
            if volume > eps:
                #hdiv_heat[hru, il] = net_heat_hru[hru] / volume
                hdiv_heat[hru, il] = net_heat_hru[hru] / area_hrus[hru]
            else:
                hdiv_heat[hru, il] = 0.0

    return hdiv_heat

@numba.jit(nopython=True,cache=True)
def calculate_soil_moisture_potential(il,theta,thetar,thetas,b,satpsi):
    eps = 0.01
    theta = theta[:,il]
    m = (theta <= (1+eps)*thetar)
    theta[m] = (1+eps)*thetar[m]
    psi = satpsi*((theta-thetar)/(thetas-thetar))**-b
    return psi

@numba.jit(nopython=True,cache=True)
def calculate_transmissivity(psi,ztop,zbot,m,ksat,satpsi,b, af):
    #af = 2.0  #Daniel
    Ksat_x = af*ksat #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
    K_x = Ksat_x*np.true_divide(psi,satpsi)**(-2-np.true_divide(3.,b))
    #Correct hydraulic conductivity if layer is below 1.5 meters
    depth_threshold = 1.5
    #Calculate transmissivity at top layer (exponential decay)
    Ttop = np.abs(np.where(zbot>depth_threshold, m*K_x*np.exp(-(ztop-depth_threshold)/m),K_x*ztop))
    #Calculate transmissivity at bottom of layer (exponential decay)
    Tbot = np.abs(np.where(zbot>depth_threshold, m*K_x*np.exp(-(zbot-depth_threshold)/m),K_x*zbot))
    # Compute transmissivity
    T = np.abs(Ttop - Tbot)
    return T

@numba.jit(nopython=True,cache=True)
def calculate_hydraulic_head(hand,psi,depth):
    h = hand - depth - psi
    return h

@numba.jit(nopython=True,cache=True)
def calculate_lateralflow_gw(h,T,w,dx,mconx):
    #Calculate dh
    dh = calculate_dh_gw(h, mconx)
    #Calculate That
    That = calculate_That_gw(T, mconx);
    #q = -That*dh/dx*w; #m3/s -1000.0*That*dh/dx*w/area (original) Daniel
    #q = That*dh/dx*w; #no need of (-) because the div is already aligned with the direction of the flow and the sum is done over axis = 1
    # Compute lateral flow element-wise to avoid division by zero when dx[i,j] = 0
    eps = 1e-20
    q = np.zeros_like(That)
    for i in range(q.shape[0]):
        for j in range(q.shape[1]):
            if dx[i,j] > eps:
                q[i,j] = That[i,j] * dh[i,j] / dx[i,j] * w[i,j]
            else:
                q[i,j] = 0.0
    return q

@numba.jit(nopython=True,cache=True)
def calculate_dh_gw(h, mconx):
    dh = np.zeros((mconx.shape))
    for i in range(h.size):
        for j in range(mconx.shape[1]):
            dh[i,j] = h[i] - h[mconx[i,j]]
    return dh

@numba.jit(nopython=True,cache=True)
def calculate_That_gw(T,mconx):
    That = np.zeros(mconx.shape)
    eps = 1e-20  # Small epsilon to avoid division by zero
    for i in range(T.size):
        for j in range(mconx.shape[1]):
            denom = T[i] + T[mconx[i,j]]
            if denom > eps:
                That[i,j] = (2*T[i]*T[mconx[i,j]]) / denom
            else:
                That[i,j] = 0.0  # If both transmissivities are ~0, set to 0
    return That
        
@numba.njit(nopython=True,cache=True)                                              
def parameters_aggregation(nl, theta, area, units):
    theta_gw = np.zeros((units,nl));
    vector2  = np.ones(area.shape[1]); #array with number of hrus
    for il in range(nl):       
        #Aggregate values to units
        vector1 = theta[:,il]*area;
        theta_gw[:,il] = sum_product(vector1,vector2)
    return theta_gw

@numba.njit(nopython=True,cache=True)                                              
def sum_product(vector1, vector2):
    result = np.zeros(vector1.shape[:-1])
    # Iterate over the last axis
    for i in range(vector1.shape[-1]):
        result += vector1[..., i] * vector2[..., i]
    return result

#-------------------Defining the functions to transfer SMC/T/DZ between regional units---------------------------    
def risfu_connections_regional(cids,rank,size,HBdb):
    # Create the mapping of connections between regional units across ranks. 
    # Each rank sends the list of rows of each cid it needs from the others to rank 0, 
    # which then builds the complete mapping and sends it back to all ranks. 
    # The mapping is stored in local_subsurface.risfu_mapping for each cid in the rank.
    cids_core = cids
    #cid_rank_mapping = HBdb[cids[0]].cid_rank_mapping
    if rank != 0:
        dest = 0
        db_ex = {}
        for cid_in_rank in cids_core:
            local_subsurface = HBdb[cid_in_rank].mssubsurface
            db_ex[cid_in_rank] = {}
            for cid in local_subsurface.reg_ids.keys():
                rows = local_subsurface.reg_ids[cid] - 1
                db_ex[cid_in_rank][cid] = rows
        local_subsurface.comm.send(db_ex,dest=dest,tag=11) #send the rows of each cid needs from the others to rank 0
    elif rank == 0:
        db = {}
        for cid_in_rank in cids_core:
            local_subsurface = HBdb[cid_in_rank].mssubsurface
            db[cid_in_rank] = {}
            for cid in local_subsurface.reg_ids.keys():
                rows = local_subsurface.reg_ids[cid] - 1
                db[cid_in_rank][cid] = rows 
                
        for i in range(1,size):
            db_ex = local_subsurface.comm.recv(source=i,tag=11)
            for key in db_ex:
                db[key] = db_ex[key]
    #Wait until completed
    local_subsurface.comm.Barrier()
    #Send the list to all the cores
    if rank == 0:
        for i in range(1,size):
            local_subsurface.comm.send(db,dest=i,tag=11)
    
    if rank != 0:
        db = local_subsurface.comm.recv(source=0,tag=11)

    #Memorize links
    for cid_in_rank in cids:
        local_subsurface = HBdb[cid_in_rank].mssubsurface
        local_subsurface.risfu_mapping = db
        #print(f'in rank {rank}: for cids: {cids} cid_in_rank: {cid_in_rank} connections: {local_subsurface.risfu_mapping}',flush=True)
    return

def exchange_smc_regional_units(cids,rank,HBdb):
    # This function performs the exchange of soil moisture content (SMC) between regional units across ranks.
    # It uses non-blocking MPI communication to send and receive the necessary SMC data based on the connections defined in the risfu_mapping. 
    # The function also includes a check to ensure that the planned sends and receives match across all ranks.
    # After the exchange, it updates the reg_theta_gw arrays with the received SMC values for the regional units.
    if len(cids) == 0:
        return

    cid_rank_mapping = HBdb[cids[0]].cid_rank_mapping
    cids_core = set(cids)
    comm = HBdb[cids[0]].mssubsurface.comm

    # Keep SMC exchanges in an exclusive tag space to avoid cross-variable collisions.
    smc_tag_offset = 200000

    request_send = []
    request_recv = []
    smc_db_all = {}

    # Track communication plans so mismatches can fail fast instead of hanging in Waitall.
    planned_send = []
    planned_recv = []

    for cid_in_rank in cids:
        local_subsurface = HBdb[cid_in_rank].mssubsurface
        smc_hb = local_subsurface.th_gw[:]
        risfu_mapping = local_subsurface.risfu_mapping

        smc_db = {cid_in_rank: smc_hb}

        for cid in risfu_mapping.keys():
            if cid not in cids_core:
                if cid_in_rank in risfu_mapping[cid].keys():
                    rows_to_send = risfu_mapping[cid][cid_in_rank]
                    data_to_send = smc_hb[rows_to_send, :]
                    dest = cid_rank_mapping[cid]
                    tag = smc_tag_offset + (cid_in_rank * 1000 + cid)

                    request_send.append(local_subsurface.comm.Isend(data_to_send, dest=dest, tag=tag))
                    planned_send.append((rank, dest, tag, data_to_send.shape[0], data_to_send.shape[1], str(data_to_send.dtype)))

            elif cid != cid_in_rank and cid in risfu_mapping[cid_in_rank].keys():
                rows = risfu_mapping[cid_in_rank][cid]
                smc_other_cid = HBdb[cid].mssubsurface.th_gw[:]
                smc_db[cid] = smc_other_cid[rows, :]

        for cid in risfu_mapping.keys():
            if cid not in cids_core and cid in risfu_mapping[cid_in_rank].keys():
                src = cid_rank_mapping[cid]
                tag = smc_tag_offset + (cid * 1000 + cid_in_rank)
                rows = risfu_mapping[cid_in_rank][cid]
                # Match dtype with local SMC arrays to avoid MPI datatype inconsistencies.
                data_recv_buffer = np.empty((len(rows), smc_hb.shape[1]), dtype=smc_hb.dtype)

                request_recv.append(local_subsurface.comm.Irecv(data_recv_buffer, source=src, tag=tag))
                planned_recv.append((src, rank, tag, data_recv_buffer.shape[0], data_recv_buffer.shape[1], str(data_recv_buffer.dtype)))
                smc_db[cid] = data_recv_buffer

        smc_db_all[cid_in_rank] = smc_db

    gathered_sends = comm.allgather(planned_send)
    gathered_recvs = comm.allgather(planned_recv)
    global_send_plan = sum(gathered_sends, [])
    global_recv_plan = sum(gathered_recvs, [])
    if Counter(global_send_plan) != Counter(global_recv_plan):
        missing_recvs = list((Counter(global_send_plan) - Counter(global_recv_plan)).elements())[:5]
        missing_sends = list((Counter(global_recv_plan) - Counter(global_send_plan)).elements())[:5]
        raise ValueError(
            f"SMC MPI plan mismatch on rank={rank}. "
            f"Missing recv matches (sample): {missing_recvs}; "
            f"missing send matches (sample): {missing_sends}"
        )

    comm.Barrier()
    all_requests = request_recv + request_send
    if len(all_requests) > 0:
        MPI.Request.Waitall(all_requests)
    comm.Barrier()

    for cid_in_rank in cids:
        local_subsurface = HBdb[cid_in_rank].mssubsurface
        smc_db = smc_db_all[cid_in_rank]

        smc = np.empty((1, local_subsurface.nsoil))
        for cid in local_subsurface.reg_ids.keys():
            if cid not in smc_db:
                raise ValueError(f"Missing smc data for cid_in_rank={cid_in_rank}, cid={cid} on rank={rank}")
            smc = np.concatenate((smc, smc_db[cid]), axis=0)

        local_subsurface.reg_theta_gw[:, :] = smc[1:, :]

    return

def exchange_temperature_regional_units(cids,rank,HBdb):
    # This function performs the exchange of groundwater temperature between regional units across ranks.
    # It uses non-blocking MPI communication to send and receive the necessary temperature data based on
    # the connections defined in the risfu_mapping.
    cid_rank_mapping = HBdb[cids[0]].cid_rank_mapping
    cids_core = cids

    request_send = []
    request_recv = []

    temperature_db_all = {}

    for cid_in_rank in cids_core:
        local_subsurface = HBdb[cid_in_rank].mssubsurface
        temperature_gw = local_subsurface.temp_gw[:]
        risfu_mapping = local_subsurface.risfu_mapping

        temperature_db = {}
        temperature_db[cid_in_rank] = temperature_gw

        for cid in risfu_mapping.keys():
            if cid not in cids_core:
                if cid_in_rank in risfu_mapping[cid].keys():
                    rows_to_send = risfu_mapping[cid][cid_in_rank]
                    data_to_send = temperature_gw[rows_to_send, :]
                    dest = cid_rank_mapping[cid]
                    tag = (cid_in_rank * 1000 + cid + 400000)
                    request = local_subsurface.comm.Isend(data_to_send, dest=dest, tag=tag)
                    request_send.append(request)

            elif cid in cids_core and cid != cid_in_rank:
                if cid in risfu_mapping[cid_in_rank].keys():
                    rows = risfu_mapping[cid_in_rank][cid]
                    temperature_other_cid = HBdb[cid].mssubsurface.temp_gw[:]
                    data_rows = temperature_other_cid[rows, :]
                    temperature_db[cid] = data_rows

        for cid in risfu_mapping.keys():
            if cid not in cids_core:
                if cid in risfu_mapping[cid_in_rank].keys():
                    src = cid_rank_mapping[cid]
                    tag = (cid * 1000 + cid_in_rank + 400000)
                    rows = risfu_mapping[cid_in_rank][cid]
                    data_recv_buffer = np.empty((len(rows), temperature_gw.shape[1]))
                    request_r = local_subsurface.comm.Irecv(data_recv_buffer, source=src, tag=tag)
                    request_recv.append(request_r)
                    temperature_db[cid] = data_recv_buffer

        temperature_db_all[cid_in_rank] = temperature_db

    HBdb[cids_core[0]].mssubsurface.comm.Barrier()
    MPI.Request.Waitall(request_send)
    MPI.Request.Waitall(request_recv)
    HBdb[cids_core[0]].mssubsurface.comm.Barrier()

    for cid_in_rank in cids_core:
        local_subsurface = HBdb[cid_in_rank].mssubsurface
        temperature_db = temperature_db_all[cid_in_rank]

        temperature = np.empty((1, local_subsurface.nsoil))
        for cid in local_subsurface.reg_ids.keys():
            if cid not in temperature_db:
                raise ValueError(f"Missing temperature data for cid_in_rank={cid_in_rank}, cid={cid} on rank={rank}")
            data_units = temperature_db[cid]
            temperature = np.concatenate((temperature, data_units), axis=0)

        local_subsurface.reg_temperature_gw[:, :] = temperature[1:, :]

    return

def exchange_dz_regional_units(cids, rank, HBdb):
    # This function performs the exchange of groundwater layer thickness (dz) between regional units across ranks.
    # It uses non-blocking MPI communication to send and receive the necessary dz data based on the connections defined in the risfu_mapping.
    if len(cids) == 0:
        return

    cid_rank_mapping = HBdb[cids[0]].cid_rank_mapping
    cids_core = set(cids)
    comm = HBdb[cids[0]].mssubsurface.comm

    # Keep dz exchanges in an exclusive tag space to avoid cross-variable collisions.
    dz_tag_offset = 300000

    request_send = []
    request_recv = []
    dz_db_all = {}

    # Track communication plans so mismatches can fail fast instead of hanging in Waitall.
    planned_send = []
    planned_recv = []

    for cid_in_rank in cids:
        local_subsurface = HBdb[cid_in_rank].mssubsurface
        dz_hb = local_subsurface.dz_gw[:]
        risfu_mapping = local_subsurface.risfu_mapping

        dz_db = {cid_in_rank: dz_hb}

        for cid in risfu_mapping.keys():
            if cid not in cids_core:
                if cid_in_rank in risfu_mapping[cid].keys():
                    rows_to_send = risfu_mapping[cid][cid_in_rank]
                    data_to_send = dz_hb[rows_to_send, :]
                    dest = cid_rank_mapping[cid]
                    tag = dz_tag_offset + (cid_in_rank * 1000 + cid)

                    request_send.append(local_subsurface.comm.Isend(data_to_send, dest=dest, tag=tag))
                    planned_send.append((rank, dest, tag, data_to_send.shape[0], data_to_send.shape[1], str(data_to_send.dtype)))

            elif cid != cid_in_rank and cid in risfu_mapping[cid_in_rank].keys():
                rows = risfu_mapping[cid_in_rank][cid]
                dz_other_cid = HBdb[cid].mssubsurface.dz_gw[:]
                dz_db[cid] = dz_other_cid[rows, :]

        for cid in risfu_mapping.keys():
            if cid not in cids_core and cid in risfu_mapping[cid_in_rank].keys():
                src = cid_rank_mapping[cid]
                tag = dz_tag_offset + (cid * 1000 + cid_in_rank)
                rows = risfu_mapping[cid_in_rank][cid]
                # Match dtype with local dz arrays to avoid MPI datatype inconsistencies.
                data_recv_buffer = np.empty((len(rows), dz_hb.shape[1]), dtype=dz_hb.dtype)

                request_recv.append(local_subsurface.comm.Irecv(data_recv_buffer, source=src, tag=tag))
                planned_recv.append((src, rank, tag, data_recv_buffer.shape[0], data_recv_buffer.shape[1], str(data_recv_buffer.dtype)))
                dz_db[cid] = data_recv_buffer

        dz_db_all[cid_in_rank] = dz_db

    gathered_sends = comm.allgather(planned_send)
    gathered_recvs = comm.allgather(planned_recv)
    global_send_plan = sum(gathered_sends, [])
    global_recv_plan = sum(gathered_recvs, [])
    if Counter(global_send_plan) != Counter(global_recv_plan):
        missing_recvs = list((Counter(global_send_plan) - Counter(global_recv_plan)).elements())[:5]
        missing_sends = list((Counter(global_recv_plan) - Counter(global_send_plan)).elements())[:5]
        raise ValueError(
            f"DZ MPI plan mismatch on rank={rank}. "
            f"Missing recv matches (sample): {missing_recvs}; "
            f"missing send matches (sample): {missing_sends}"
        )

    comm.Barrier()
    all_requests = request_recv + request_send
    if len(all_requests) > 0:
        MPI.Request.Waitall(all_requests)
    comm.Barrier()

    for cid_in_rank in cids:
        local_subsurface = HBdb[cid_in_rank].mssubsurface
        dz_db = dz_db_all[cid_in_rank]

        dz = np.empty((1, local_subsurface.nsoil))
        for cid in local_subsurface.reg_ids.keys():
            if cid not in dz_db:
                raise ValueError(f"Missing dz data for cid_in_rank={cid_in_rank}, cid={cid} on rank={rank}")
            dz = np.concatenate((dz, dz_db[cid]), axis=0)

        local_subsurface.reg_dz_gw[:, :] = dz[1:, :]

    return