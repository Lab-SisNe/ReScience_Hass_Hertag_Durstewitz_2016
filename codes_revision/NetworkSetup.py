from BasicsSetup import basics_setup, MembraneTuple, BaseClass
import numpy as np
import xarray as xr
from dataclasses import dataclass
from scipy.integrate import quad
from scipy.optimize import fmin, fsolve
from warnings import filterwarnings
import os, shutil, json
from importlib import import_module
from Auxiliar import time_report

filterwarnings("ignore", category=RuntimeWarning) 

@time_report('Network setup')
def network_setup(Ncells_prompted, Nstripes, basics_scales=None, seed=None):   
    
    np.random.seed(seed)

    basics=basics_setup(Ncells_prompted, Nstripes, basics_scales)
    
    index_stripes = np.arange(Nstripes).astype(int)
    index_cells_total = np.arange(basics.struct.Ncells_total).astype(int)
    
    syn_pairs = xr.DataArray(np.asarray([[],[]]), coords=[['target', 'source'],[]], dims=['cell', 'syn_index'])
    
    syn_channel = ['channel']
    syn_spike_params = ['gmax', 'pfail']
    syn_STSP_params = ['U', 'tau_rec', 'tau_fac']
    syn_delay = ['delay']
    
    syn_params = {}
    syn_params['channel'] = np.asarray([[] for i in range(len(syn_channel))])
    syn_params['spiking'] =np.asarray([[] for i in range(len(syn_spike_params))])
    syn_params['STSP_params'] = np.asarray([[] for i in range(len(syn_STSP_params))])
    syn_params['delay'] = np.asarray([[] for i in range(len(syn_delay))])
    # synparams = xr.DataArray(np.asarray([[] for i in range(len(synparams))]), coords=[synparams, []], dims=['param', 'syn_index'])
    
    membr_params = xr.DataArray(np.zeros((len(basics.membr.names), basics.struct.Ncells_total)), coords=(basics.membr.names, index_cells_total), dims=['param', 'cell_index'])
    I_refractory = xr.DataArray(np.zeros(len(index_cells_total)), coords=[index_cells_total,], dims='cell_index')
    
    # -------------- Set neuron parameters and most synaptic connections -----------------
 
    group_distr = [ [[] for j in range(basics.struct.groups.N)] for i in range (Nstripes)]

    Nsyn_current = 0
    Ncell_current = 0
    for stripe in index_stripes:   
        for group in basics.struct.groups.names:
            Ncell_new=int(basics.struct.Ncells_per_group.loc[group].values)
            set_current = np.arange(Ncell_new).astype(int)+Ncell_current
            Ncell_current+=Ncell_new
          
            group_distr[stripe][basics.struct.groups.names.index(group)].extend(set_current)
            set_new = set_current
            
            for trial in range(1000):
                multi_mean = np.zeros(len(basics.membr.names))
                multi_cov = basics.membr.covariance.loc[dict(group=group)]
                multi_N = len(set_new)
                
                NeuPar_multi = np.random.multivariate_normal(multi_mean,multi_cov,multi_N).transpose()
                NeuPar_multi=xr.DataArray(NeuPar_multi, coords=[basics.membr.names, set_new], dims=['param', 'cell_index'])
                
                membr_params, set_new = set_membr_params(membr_params, set_new, set_current, 
                                                          NeuPar_multi, group, basics)  

                if len(set_new)==0:                             
                    break

            else:
                print('WARNING: Number of trials for the full multivariate distribution has been exceeded. Use multivariate uniform distribution for the rest.\n')
                
                for trial in range(1000):
                    
                    NeuPar_unif=np.random.uniform(basics.membr.min.loc[dict(group=group)],basics.membr.max.loc[dict(group=group)], len(set_new))
                    NeuPar_unif=xr.DataArray(NeuPar_unif, coords=[basics.membr.names, set_new], dims=['param', 'cell_index'])

                    membr_params, set_new = set_membr_params(membr_params, NeuPar_unif, group, basics)
                    
                    if len(set_new)==0:
                        break
                
                else:
                    print('ERROR: Number of trials for univariate uniform distribution has been exceeded.\n')
                         

            for cell in set_current:#index_cells_total:
                memb = MembraneTuple(*membr_params.loc[dict(cell_index=cell)].to_numpy())                
                Irheo = memb.g_L * (memb.V_T - memb.E_L) - memb.g_L * memb.delta_T                 
                I_refractory[cell] = fmin(Define_I_ref, Irheo+100, disp=False, args=(memb,))[0]
                
        group_distr = redistribute(basics, group_distr, membr_params, stripe)
         
        #     # intra-celltype connections
        
        for group in basics.struct.groups.names:
            set_current = group_distr[stripe][basics.struct.groups.names.index(group)]
           
            pCon = float(basics.struct.conn.pCon.loc[group, group].values)
            if len(set_current)>0:
                if basics.struct.conn.cluster.loc[group, group]==1:
                    TS_pairs, syn_idc = SetCon_CommonNeighbour(Nsyn_current, set_current, pCon, 0.47)
                    
                else: 
                    TS_pairs, syn_idc = SetCon_standard(Nsyn_current, set_current, set_current, pCon)         # ... or without common neighbour rule
            
                if len(syn_idc)>0:
                    syn_params, syn_pairs, Nsyn_current = set_syn_params(basics, 1, 1, group, group, 
                                                                         syn_idc, TS_pairs, syn_params, syn_pairs, Nsyn_current)
                    
                
        for group_target in basics.struct.groups.names:
            set_current_target=group_distr[stripe][basics.struct.groups.names.index(group_target)]
            
            for group_source in [name for name in basics.struct.groups.names if name!=group_target]:
                set_current_source=group_distr[stripe][basics.struct.groups.names.index(group_source)]
                pCon = float(basics.struct.conn.pCon.loc[group_target, group_source].values)
                
                if len(set_current_target)*len(set_current_source)>0:
        
                    TS_pairs, syn_idc = SetCon_standard(Nsyn_current, set_current_target, set_current_source,pCon)         # ... or without common neighbour rule
    
                    if TS_pairs.shape[1]>0:
                        syn_params, syn_pairs, Nsyn_current = set_syn_params(basics, 1, 1, group_target, group_source, 
                                                                             syn_idc, TS_pairs, syn_params, syn_pairs, Nsyn_current)
                        
    
    # # # ------------------------- Define inter-stripe connections ---------------------------------
    
    if Nstripes>1:
        
        for intstp_set in basics.struct.stripes.inter:
            intstp_dict = basics.struct.stripes.inter[intstp_set]
            
            group_target, group_source  =intstp_dict['pair']
            pCon  = float(basics.struct.conn.pCon.loc[group_target, group_source].values)
            
            
            for stripe in range(Nstripes):             
                for conn  in range(len(intstp_dict['connections'])):
                    target_stripe = stripe + intstp_dict['connections'][conn]
                    
                    while target_stripe<0:
                        target_stripe += Nstripes
                    
                    while target_stripe>=Nstripes:
                        target_stripe -= Nstripes
                    
                    dist_act = abs(intstp_dict['connections'][conn])
                    pCon_curr = pCon*np.exp(-dist_act/intstp_dict['coefficient_0'])
                    Ntarget = len(group_distr[target_stripe][basics.struct.groups.names.index(group_target)]) 
                    Nsource = len(group_distr[target_stripe][basics.struct.groups.names.index(group_source)]) 
                    
                    TS_pairs,syn_idc = SetCon_standard(Nsyn_current, Ntarget, Nsource, pCon_curr)                                                                                                     # ... or without common neighbour rule
                    
                    if len(syn_idc)>0:
                        gmax_fac = np.exp(-dist_act/intstp_dict['coefficient_0'])
                        delay_fac = intstp_dict['coefficient_1']*dist_act
                                          
                        syn_params, syn_pairs, Nsyn_current = set_syn_params(basics, gmax_fac, delay_fac, group_target, group_source, 
                                                                             syn_idc, TS_pairs, syn_params, syn_pairs, Nsyn_current)
                        
    syn_params['channel'] = xr.DataArray(syn_params['channel'], coords=[syn_channel, np.arange(syn_pairs.shape[1])], dims=['param', 'syn_index'])
    syn_params['spiking'] = xr.DataArray(syn_params['spiking'], coords=[syn_spike_params,  np.arange(syn_pairs.shape[1])], dims=['param', 'syn_index'])
    syn_params['STSP_params'] = xr.DataArray(syn_params['STSP_params'], coords=[syn_STSP_params,  np.arange(syn_pairs.shape[1])], dims=['param', 'syn_index'])
    syn_params['delay'] = xr.DataArray(syn_params['delay'], coords=[syn_delay,  np.arange(syn_pairs.shape[1])], dims=['param', 'syn_index'])
    #
                        
    network = Network(membr_params, I_refractory, syn_params, syn_pairs.astype(int), group_distr, seed, basics)
    
    return  network

def set_syn_params(basics, gmax_fac, delay_fac, group_target, group_source, 
                   syn_idc, TS_pairs, syn_params, syn_pairs, Nsyn_current):
    syntypes_current = str(basics.syn.kinds.loc[group_target, group_source].values)
    synchannels_current = basics.syn.channels.kinds_to_names[syntypes_current]
     
    gmax = float(basics.syn.spiking.params.gmax['mean'].loc[group_target, group_source])*gmax_fac
    gmax_sigma = float(basics.syn.spiking.params.gmax['sigma'].loc[group_target, group_source])*gmax_fac
    gmax_min = float(basics.syn.spiking.params.gmax['min'].loc[group_target, group_source])
    gmax_max = float(basics.syn.spiking.params.gmax['max'].loc[group_target, group_source])
    
    delay = float(basics.syn.delay.delay_mean.loc[group_target, group_source])+delay_fac
    delay_sigma = float(basics.syn.delay.delay_sigma.loc[group_target, group_source])+delay_fac
    delay_min = float(basics.syn.delay.delay_min.loc[group_target, group_source])
    delay_max = float(basics.syn.delay.delay_max.loc[group_target, group_source])
    
    STSPset =str(basics.syn.STSP.distr.loc[group_target, group_source].values)          
    
    STSPwhere = np.isfinite((np.array(basics.syn.STSP.sets.loc[dict(set=STSPset)].values)))
    STSPkinds = np.asarray(basics.syn.STSP.sets.coords['kind'].values)[STSPwhere]
    STSPvalues = np.asarray(basics.syn.STSP.sets.loc[dict(set=STSPset)].values)[STSPwhere]

    new_delay=set_syn_delay(syn_idc, delay, delay_sigma, delay_min, delay_max)

    for channel in synchannels_current:
        syn_idc = syn_idc.astype(int)
        new_channel = xr.DataArray([[basics.syn.channels.names.index(channel)]*len(syn_idc)], coords=[['channel',], syn_idc], dims=['param', 'syn_index'])      
        new_STSP_params = set_syn_STSP(syn_idc, STSPkinds,STSPvalues, basics.syn.STSP.kinds)
        new_gmax = set_syn_gmax(syn_idc, gmax, gmax_sigma, gmax_min, gmax_max)
        new_pfail = set_syn_pfail(syn_idc,  basics.syn.spiking.params.pfail.loc[group_target, group_source])
        new_spike_params = np.concatenate((new_gmax, new_pfail), axis=0)
        new_delay=new_delay.assign_coords(syn_index=syn_idc)
        
        syn_params['channel'] = np.concatenate((syn_params['channel'], new_channel), axis=1)
        syn_params['spiking'] = np.concatenate((syn_params['spiking'], new_spike_params), axis=1)
        syn_params['STSP_params'] = np.concatenate((syn_params['STSP_params'], new_STSP_params), axis=1)
        syn_params['delay'] = np.concatenate((syn_params['delay'], new_delay), axis=1)
        
        syn_pairs = np.concatenate((syn_pairs, TS_pairs), axis=1)
        Nsyn_current +=len(syn_idc)        

    return syn_params,  syn_pairs, Nsyn_current, 

def redistribute(basics, group_distr, membr_params, stripe):
    
    I0 = 500

    set_current = []
    for group in ['IN_L_L23', 'IN_L_d_L23']:         
            set_current.extend(group_distr[stripe][basics.struct.groups.names.index(group)])

    set_L = []
    set_L_d = []
    for cell in set_current:
        memb = MembraneTuple(*membr_params[dict(cell_index=cell)].to_numpy())
        t_lat_AdEx = latency_AdEx(memb, I0)
        t_lat_LIF =  latency_LIF(memb, I0)
        if t_lat_AdEx-t_lat_LIF>0:
            set_L_d.append(cell)
        else:
            set_L.append(cell)

    group_distr[stripe][basics.struct.groups.names.index('IN_L_L23')] = set_L
    group_distr[stripe][basics.struct.groups.names.index('IN_L_d_L23')] = set_L_d
            
    set_current = []
    for group in ['IN_L_L5', 'IN_L_d_L5']:         
            set_current.extend(group_distr[stripe][basics.struct.groups.names.index(group)])

    set_L = []
    set_L_d = []
    for cell in set_current:
        memb = MembraneTuple(*membr_params.loc[dict(cell_index=cell)].to_numpy())         
        t_lat_AdEx = latency_AdEx(memb, I0)
        t_lat_LIF =  latency_LIF(memb, I0)
        if t_lat_AdEx-t_lat_LIF>0:
            set_L_d.append(cell)
        else:
            set_L.append(cell)

    group_distr[stripe][basics.struct.groups.names.index('IN_L_L5')] = set_L
    group_distr[stripe][basics.struct.groups.names.index('IN_L_d_L5')] = set_L_d
    
    I_range = np.arange(25,301,25)
    adapt_cut = 1.5834
    
    set_current = []    
    for group in ['IN_CL_L23', 'IN_CL_AC_L23']:
        set_current.extend(group_distr[stripe][basics.struct.groups.names.index(group)])
    
    set_CL = []
    set_CL_AC = []
    for cell in set_current:  
        f_trans = np.zeros(len(I_range))
        f_stat  = np.zeros(len(I_range))            
        for Iidc in range(len(I_range)):
            memb_par = MembraneTuple(*membr_params.loc[dict(cell_index=cell)].to_numpy())
            f_trans[Iidc] = get_transFR(memb_par,I_range[Iidc],0,memb_par.E_L)
            f_stat[Iidc]  = get_statFR(memb_par,I_range[Iidc])
                    
        idc_valid = np.where((f_stat>0))[0]  
        if np.median(f_trans[idc_valid] / f_stat[idc_valid])>adapt_cut:
            set_CL_AC.append(cell)
        else:
            set_CL.append(cell)
      
    group_distr[stripe][basics.struct.groups.names.index('IN_CL_L23')] = set_CL
    group_distr[stripe][basics.struct.groups.names.index('IN_CL_AC_L23')] = set_CL_AC
    
    set_current = []    
    for group in ['IN_CL_L5', 'IN_CL_AC_L5']:
        set_current.extend(group_distr[stripe][basics.struct.groups.names.index(group)])
    
    set_CL = []
    set_CL_AC = []
    for cell in set_current:  
        f_trans = np.zeros(len(I_range))
        f_stat  = np.zeros(len(I_range))            
        for Iidc in range(len(I_range)):
            memb_par = MembraneTuple(*membr_params.loc[dict(cell_index=cell)].to_numpy())
            f_trans[Iidc] = get_transFR(memb_par,I_range[Iidc],0,memb_par.E_L)
            f_stat[Iidc]  = get_statFR(memb_par,I_range[Iidc])
                    
        idc_valid = np.where((f_stat>0))[0]  
        if np.median(f_trans[idc_valid] / f_stat[idc_valid])>adapt_cut:
            set_CL_AC.append(cell)
        else:
            set_CL.append(cell)
      
    group_distr[stripe][basics.struct.groups.names.index('IN_CL_L5')] = set_CL
    group_distr[stripe][basics.struct.groups.names.index('IN_CL_AC_L5')] = set_CL_AC
    
    return group_distr

def set_membr_params(membr_params, set_new, set_current, membr_par_distr, group, basics):
    
    for param in basics.membr.names:
        multi_curr = membr_par_distr.loc[dict(param=param)]
        k_curr = basics.membr.k.loc[dict(param=param, group=group)]
        mean_curr = basics.membr.mean.loc[dict(param=param, group=group)]
        std_curr = basics.membr.std.loc[dict(param=param, group=group)]
        min_curr = basics.membr.min.loc[dict(param=param, group=group)]
        inv_curr = inv_transform(multi_curr,k_curr, mean_curr, std_curr, min_curr)                  
        membr_params.loc[dict(param=param, cell_index=set_new)] = inv_curr

    membr_params.loc['C',set_new] =  membr_params.loc['C',set_new] * membr_params.loc['g_L',set_new]    # C from g_L and tau
    set_outer = np.asarray([])
    
    for param in basics.membr.names:
        where_minmax = np.where((membr_params.loc[dict(param=param,cell_index=set_current)]<basics.membr.min.loc[dict(param=param,group=group)]) | (membr_params.loc[dict(param=param,cell_index=set_current)] > basics.membr.max.loc[dict(dict(param=param,group=group))]))[0]
        set_outer = np.concatenate((set_outer, where_minmax))
    
    where_V = np.where((membr_params.loc['V_r',set_current]>=membr_params.loc['V_T',set_current]))[0]   # Vr must be smaller than Vth
    
    where_taum = np.where((membr_params.loc['C',set_current]/membr_params.loc['g_L',set_current]<basics.membr.tau_m_min.loc[group]) |(membr_params.loc['C',set_current]/membr_params.loc['g_L',set_current]>basics.membr.tau_m_max.loc[group]))[0] # check tau boundaries)
    where_tauw = np.where((membr_params.loc['tau_w',set_current]<=membr_params.loc['C',set_current]/membr_params.loc['g_L',set_current]))[0]
    where_nan = np.argwhere(np.sum(np.isnan(membr_params.loc[dict(cell_index=set_current)]), axis=0).to_numpy())[:,0]
   
    set_outer = np.concatenate((set_outer,where_V, where_taum, where_tauw, where_nan))
    set_new = set_current[np.unique(set_outer).astype(int)]  

    return membr_params, set_new

def SetCon(Nsyn, target_arr, source_arr, pCon):

    
    Ntarget, Nsource = len(target_arr), len(source_arr)
    conn_idc = np.arange(Ntarget*Nsource)
    np.random.shuffle(conn_idc)
    conn_idc=conn_idc[:int(np.round(pCon*Ntarget*Nsource,0))]
    syn_idc = np.arange(int(np.round(pCon*Ntarget*Nsource,0)))+Nsyn
   
    return get_TSpairs(conn_idc, Nsource), syn_idc.astype(int)

def SetCon_standard(Nsyn, target_arr, source_arr, pCon): 
    target_arr = np.asarray(target_arr)
    source_arr = np.asarray(source_arr)
    TS_pairs, syn_idc = SetCon(Nsyn, target_arr, source_arr, pCon)   
    return np.array([target_arr[TS_pairs[0,:]], source_arr[TS_pairs[1,:]]]), syn_idc.astype(int)


def SetCon_CommonNeighbour(Nsyn, cells_arr, pCon, pSelfCon):
    
    def get_tril(a, output='pairs'):
        if output=='pairs':
            return a[:, np.where((a[0,:]>a[1,:]))[0]]
        elif output=='indices':
            return np.where((a[0,:]>a[1,:]))[0]


    def get_bidirect(a):
        temp_a = a.copy()
        inverted = a.copy()
        inverted[0,:], inverted[1,:] = temp_a[1,:],temp_a[0,:]
        concatenated = np.concatenate((a,inverted),axis=1)
        return np.unique(concatenated, axis=1)


    def get_connections(TS_pairs, Ncells):
        connections = {}
        
        X = np.zeros((Ncells, Ncells))
        X[TS_pairs[0,:], TS_pairs[1,:]] = 1
        Xnondiag0 = (np.ones((Ncells,Ncells))-np.identity(Ncells))*X
        Xnondiag1 = X + np.identity(Ncells)
        
        Xrecur0 = Xnondiag0 + Xnondiag0.transpose()
        Xrecur1 = Xnondiag1 + Xnondiag1.transpose()
        
        connections['pairs_connected_recur'] = np.array(np.where(Xrecur0==2))
        connections['pairs_connected_all'] = np.array(np.where(Xrecur0>=1))
        connections['pairs_disconnected'] = np.array(np.where(Xrecur1==0))
        Xnondiag0[connections['pairs_connected_recur'][0,:], connections['pairs_connected_recur'][1,:]]=0
        connections['pairs_connected_uni'] = np.array(np.where(Xnondiag0==1))
       
        Xdiag0 = X * np.identity(Ncells)
        Xdiag1 = X + (np.ones((Ncells,Ncells))-np.identity(Ncells))
        connections['diag_connected'] = np.array(np.where(Xdiag0==1))
        connections['diag_disconnected'] = np.array(np.where(Xdiag1==0))
        
        return connections 


    def get_intersect(a, b, Ncells):
        a_mat = np.zeros((Ncells, Ncells))
        b_mat = a_mat.copy()
        a_mat[a[0,:],a[1,:]]=1
        b_mat[b[0,:],b[1,:]]=1
        return np.array(np.where(a_mat*b_mat>0))    

    def get_diff(a,b):
        
        if a.shape[1]==0:
            return np.array([[],[]])
        elif b.shape[1]==0:
            return a
        
        a = np.array(a).astype(int)
        b = np.array(b).astype(int)
        
        Ncells = max(np.max(a), np.max(b))+1
        a_mat = np.zeros((Ncells, Ncells))
        b_mat = a_mat.copy()
        a_mat[a[0,:],a[1,:]]=1
        b_mat[b[0,:],b[1,:]]=1
        
        return np.array(np.where(a_mat-b_mat==1))

    
    def p_calc_recur(TS_pairs, Neigh_mat, pCon, pSelfCon, slope):
        
        def p_min(N0, N_neighbours, all_neigh_count, pCon, slope):

            N1 = np.round(min(max(N_neighbours), N0 + 1/slope/pCon))
            
            N_2 = N_neighbours[np.where((N_neighbours>N0)&(N_neighbours<=N1))]
            pN_2 = all_neigh_count[np.where((N_neighbours>N0)&(N_neighbours<=N1))]/np.sum(all_neigh_count)
            pN_3 = all_neigh_count[np.where(N_neighbours>N1)]/np.sum(all_neigh_count)
            
            return abs(np.sum(pN_2 * slope * (N_2-N0)) + np.sum(pN_3) - 1)
        
        Ncell = Neigh_mat.shape[0]    
        N_neigh_connected = Neigh_mat[TS_pairs[0,:], TS_pairs[1, :]]  
        
        N_neighbours = np.unique(Neigh_mat)
        all_neigh_count = np.zeros(len(N_neighbours))
        neigh_factor = np.zeros(len(N_neighbours))
        connected_neigh_count = np.zeros(len(N_neighbours))
        for i in range(len(N_neighbours)):
            all_neigh_count[i] = np.sum((Neigh_mat == N_neighbours[i]).astype(int))
            connected_neigh_count[i] = np.sum((N_neigh_connected==N_neighbours[i]).astype(int))

        offset = fmin(lambda N0: p_min(N0, N_neighbours, all_neigh_count, pCon, slope), np.max(N_neighbours), disp=False)
         
        N1 = np.floor(min(np.max(N_neighbours), offset+1/slope/pCon ))
        N0 = max(np.min(N_neighbours), np.ceil(offset))
        ind = np.isin(N_neighbours, np.arange(N0, N1+1))
        neigh_factor[ind] = pCon*slope*(N_neighbours[ind]-offset)
        
        if len(ind)>0 and len(np.where(ind)[0])>0:
            neigh_factor[np.where(ind)[0][-1]+1:] = 1
        else:
            neigh_factor[:] = 1
        
        N_original = Ncell **2 * pCon - Ncell * pSelfCon

        N_act = all_neigh_count @ neigh_factor
        N_neigh_connections =np.round(neigh_factor*N_original*all_neigh_count/N_act).astype(int)
        
        return N_neighbours.astype(int), N_neigh_connections, all_neigh_count, connected_neigh_count
        
    cells_arr = np.array(cells_arr)
    Ncells = len(cells_arr)
    TS_pairs, syn_idc = SetCon(Nsyn, cells_arr, cells_arr, pCon)    
    Neigh_mat = get_commonneigh_recur(TS_pairs, Ncells)
    slope = 20*3.9991 / Ncells
    
    N_neighbours, N_neigh_connections, all_neigh_count, connected_neigh_count  = p_calc_recur(TS_pairs, Neigh_mat, pCon, pSelfCon, slope)

    connections= get_connections(TS_pairs, Ncells)
         
    curr_N_diag = int(round(pSelfCon*Ncells))
    
    pairs_select = np.array([[],[]])
    
    
    for n in [n for n in range(len(N_neighbours)) if N_neigh_connections[n]>0]: # para cada número de conexões
        
        curr_all_pairs = np.array(np.where(Neigh_mat==N_neighbours[n]))  # pair idc 
        curr_pairs_disconnected = get_intersect(connections['pairs_disconnected'], curr_all_pairs, Ncells)
        curr_pairs_disconnected_tril = get_tril(curr_pairs_disconnected)
        curr_pairs_connected_recur = get_intersect(connections['pairs_connected_recur'], curr_all_pairs, Ncells)
        
        curr_pairs_connected_recur_tril = get_tril(curr_pairs_connected_recur)
 
        curr_pairs_connected_uni = get_intersect(connections['pairs_connected_uni'], curr_all_pairs, Ncells)
        curr_N_recur = np.floor(N_neigh_connections[n]*pSelfCon/2).astype(int)# N recorr
        curr_N_uni = N_neigh_connections[n] - 2 * curr_N_recur # N uni
        
        if curr_pairs_connected_recur_tril.shape[1] >= curr_N_recur: # mais conexões recorrentes que o estipulado
            idc_curr= np.arange(curr_pairs_connected_recur_tril.shape[1])
            np.random.shuffle(idc_curr)        
 
            idc_select = np.asarray(idc_curr[:curr_N_recur]) # tomando N_rec aleatórios      
            pairs_select = np.concatenate((pairs_select, get_bidirect(curr_pairs_connected_recur_tril[:, idc_select])), axis=1)
            
        else: # menos conexões recorrentes que o estipulado
           
            pairs_select = np.concatenate((pairs_select, curr_pairs_connected_recur), axis=1)
                    
            N_new_idc = curr_N_recur  - curr_pairs_connected_recur_tril.shape[1] # quantas faltam
            
            idc_new = np.arange(curr_pairs_disconnected_tril.shape[1])
            np.random.shuffle(idc_new)
            idc_select = np.asarray(idc_new[:N_new_idc]) # tomando N_new_idc aleatórias k1_new
            pairs_select = np.concatenate((pairs_select, get_bidirect(curr_pairs_disconnected_tril[:, idc_select])), axis=1)
            curr_pairs_disconnected = get_diff(curr_pairs_disconnected, get_bidirect(curr_pairs_disconnected_tril[:, idc_select]))


         
        if curr_pairs_connected_uni.shape[1] > curr_N_uni:   # N uni antigos maior que o estipoulado   

            idc_curr = np.arange(curr_pairs_connected_uni.shape[1])
            np.random.shuffle(idc_curr)
            idc_select = np.asarray(idc_curr[:curr_N_uni]) # Tomar N_uni aleatórios
            pairs_select = np.concatenate((pairs_select, curr_pairs_connected_uni[:, idc_select]), axis=1) #Tomar N_uni aleatórios  
        
        else: # N uni antigos menor que o estipulado
            pairs_select = np.concatenate((pairs_select, curr_pairs_connected_uni), axis=1)
            N_new_idc = curr_N_uni - curr_pairs_connected_uni.shape[1] # quantos falatam

            idc_new = np.arange(curr_pairs_disconnected.shape[1])
            np.random.shuffle(idc_new)

            idc_select = np.asarray(idc_new[:N_new_idc])
            pairs_select = np.concatenate((pairs_select, curr_pairs_disconnected[:, idc_select]), axis=1) #tomando N_new_idc aleatórios
        
    curr_diag_connected = connections['diag_connected']
    curr_diag_disconnected = connections['diag_disconnected']
    
    if curr_diag_connected.shape[1] > curr_N_diag: # se maior que o estipulado
        idc_old = np.arange(curr_diag_connected.shape[1])
        np.random.shuffle(idc_old)
        idc_select = idc_old[:curr_N_diag]
        pairs_select = np.concatenate((pairs_select, curr_diag_connected[:, idc_select]), axis=1)   # tomando apenas o n estiulado aleatoriamente     
    else:
        N_new_idc = curr_N_diag - curr_diag_connected.shape[1] # quantos faltam
        
        pairs_select = np.concatenate((pairs_select, curr_diag_connected), axis=1)
    
        idc_new = np.arange(curr_diag_disconnected.shape[1])
        np.random.shuffle(idc_new)
        idc_select = idc_new[:N_new_idc]

        pairs_select  = np.concatenate((pairs_select, curr_diag_disconnected[:,idc_select]), axis=1) # aleatoriamente
        
    idc = Nsyn + np.arange(pairs_select.shape[1]) # gerando novo idc
    
    pairs_select = pairs_select.astype(int)

    return np.array([cells_arr[pairs_select[0,:]], cells_arr[pairs_select[1,:]]]), idc.astype(int)



def get_TSpairs(conn_idc, Nsource):
    return np.asarray([conn_idc//Nsource,conn_idc%Nsource]).astype(int)

        
def set_syn_gmax(new_idc, gmax, gmax_sigma, gmax_min, gmax_max):
       
    syn_params = xr.DataArray(np.zeros((1, len(new_idc))), coords=[['gmax'], new_idc], dims=['param', 'syn_index'])

    mean_gmax = np.log((gmax**2)/np.sqrt(gmax_sigma**2 + gmax **2))
    std_gmax = np.sqrt(np.log(gmax_sigma**2/gmax**2 + 1))
    
    syn_params.loc['gmax', new_idc] = random_param(len(new_idc), mean_gmax, std_gmax, gmax_min*gmax, gmax_max * gmax, 'log_normal')
           
    return syn_params

def set_syn_pfail(new_idc, pfail):
       
    syn_params = xr.DataArray(np.zeros((1, len(new_idc))), coords=[['pfail'], new_idc], dims=['param', 'syn_index'])
    syn_params.loc['pfail', new_idc]=pfail
    
    return syn_params

def set_syn_delay(new_idc, delay, delay_sigma, delay_min, delay_max):
   
    syn_params = xr.DataArray(np.zeros((1, len(new_idc))), coords=[['delay'], new_idc], dims=['param', 'syn_index'])
    syn_params.loc['delay', new_idc]= random_param(len(new_idc), delay, delay_sigma, delay_min*delay, delay_max*delay, 'normal')
    
    return syn_params

def set_syn_STSP(new_idc, STSPkinds, STSPvalues, STSPparams):
    
    
    syn_params = xr.DataArray(np.zeros((3, len(new_idc))), coords=[['U', 'tau_rec', 'tau_fac'], new_idc], dims=['param', 'syn_index'])
    
    # STSPkinds = np.asarray(list(STSPset.keys()))
    # STSPvalues = np.asarray(list(STSPset.values()))
    STSPidc = np.round(len(new_idc)*STSPvalues,0).astype(int)
    
    while sum(STSPidc)!=len(new_idc):   
        if sum(STSPidc) > len(new_idc):
            change_idc = np.argmin(len(new_idc)*STSPvalues - np.floor(len(new_idc)*STSPvalues))
            STSPidc[change_idc] -= 1
        else:
            change_idc = np.argmax(len(new_idc)*STSPvalues - np.floor(len(new_idc)*STSPvalues))
            STSPidc[change_idc] += 1
    
    STSPidc_cumm = np.cumsum(STSPidc).astype(int)
    STSPidc_limits = np.concatenate([np.asarray([0]), STSPidc_cumm])+new_idc[0]

    for name_idc in range(len(STSPkinds)):
        STSPpar = STSPparams.loc[dict(kind=STSPkinds[name_idc])]
        idc_curr = np.arange(STSPidc_limits[name_idc], STSPidc_limits[name_idc+1])

        syn_params.loc['U', idc_curr] = random_param(len(idc_curr), float(STSPpar.loc['U_mean'].values), float(STSPpar.loc['U_std'].values),    0,    1, 'normal')   
        syn_params.loc['tau_rec', idc_curr] = random_param(len(idc_curr), float(STSPpar.loc['tau_rec_mean'].values), float(STSPpar.loc['tau_rec_std'].values), 0, 1500, 'normal')
        syn_params.loc['tau_fac', idc_curr] = random_param(len(idc_curr), float(STSPpar.loc['tau_fac_mean'].values), float(STSPpar.loc['tau_fac_std'].values), 0, 1500, 'normal')

    return syn_params

                            
def inv_transform(X_trans, k, mean_X, std_X, min_X):
    
    std_decr =0.8
    
    if k>0:
        if min_X<0:
            X_inv = (mean_X+std_decr*std_X*X_trans)**(1/k)+1.1*min_X
        else:
            X_inv = (mean_X+std_decr*std_X*X_trans)**(1/k)
    
    else:
        if min_X<0:
            X_inv = np.exp(mean_X+std_decr*std_X*X_trans)+1.1*min_X;
        else:
            X_inv = np.exp(mean_X+std_decr*std_X*X_trans)

    return X_inv

def get_statFR(memb, I0):
    
    tau_m = memb.C / memb.g_L
    tau_ratio=tau_m / memb.tau_w
    Dist_VT = tau_ratio*(I0 + memb.g_L * memb.delta_T - memb.g_L * (memb.V_T - memb.E_L))                  
    w_end = -memb.g_L * (memb.V_T - memb.E_L) + memb.g_L * memb.delta_T + I0 - Dist_VT
  
    if memb.b!=0:
        w_r=w_end+memb.b
    else:
        w_r=0    
        
    return get_FR(memb, I0, w_r, memb.V_r)

def get_transFR(memb, I0, w0,V0):

        
    return get_FR(memb, I0, w0, V0)

def get_FR(memb, I0, w_r, V_r):
    
    def ISI_integrate_part1and3(V,I,memb, w_r):
        
        F=(1/memb.C)*(I-w_r+memb.g_L*memb.delta_T*np.exp((V-memb.V_T)/memb.delta_T)-memb.g_L*(V-memb.E_L))
        f=1/F
        
        return f

    def ISI_integrate_part2(V,I0,memb):
        
        tau_ratio=memb.C / (memb.tau_w*memb.g_L)
        
        k0 = (tau_ratio-1)*memb.g_L
        k1 = (1-tau_ratio)*I0-k0*memb.E_L
        k2=(1-tau_ratio)*memb.g_L*memb.delta_T
        
        F=(1/memb.C)*(I0-(k0*V+k1+k2*np.exp((V-memb.V_T)/memb.delta_T))+memb.g_L*memb.delta_T*np.exp((V-memb.V_T)/memb.delta_T)-memb.g_L*(V-memb.E_L))
        f=1/F
        
        return f

    def get_w_nullbound_inters(memb,I,w0,w_ref,signal):

        
        if w0>w_ref:
            V_nullbound_inters=[None, None]
            tau_ratio=memb.C/(memb.g_L*memb.tau_w)
         
            nullcl_bound= lambda V: (1+signal*tau_ratio)*(I-memb.g_L*(V-memb.E_L)+memb.g_L*memb.delta_T*np.exp((V-memb.V_T)/memb.delta_T)) - w0
            lo_guess=memb.E_L+(I-(w0/(1+signal*tau_ratio)))/memb.g_L-0.1
            hi_guess=memb.E_L+memb.delta_T+(I-(w0/(1+signal*tau_ratio)))/memb.g_L
            
            for trial in range(1000):         
                if np.sign(nullcl_bound(hi_guess))*np.sign(nullcl_bound(lo_guess))==-1:
                    break
                lo_guess-=1
            else:
                print('Error in numcor')
            
         
            V_nullbound_inters[0] =fsolve(nullcl_bound,(lo_guess+hi_guess)/2)
           
            lo_guess=memb.E_L+memb.delta_T+(I-(w0/(1+signal*tau_ratio)))/memb.g_L
            hi_guess=memb.V_up
            
            for trial in range(1000):
                if np.sign(nullcl_bound(hi_guess))*np.sign(nullcl_bound(lo_guess))==-1:
                    break
                lo_guess-=1
            else:
                print('Error')
            
            V_nullbound_inters[1]= fsolve(nullcl_bound,(lo_guess+ hi_guess)/2)    
            
        elif w0==w_ref:
            V_nullbound_inters=[memb.V_T,]      
        else:
            V_nullbound_inters=[] 
        
        return V_nullbound_inters


    tau_m = memb.C / memb.g_L
    tau_ratio=tau_m / memb.tau_w
    Dist_VT = tau_ratio*(I0 + memb.g_L * memb.delta_T - memb.g_L * (memb.V_T - memb.E_L))                  
    w_end = -memb.g_L * (memb.V_T - memb.E_L) + memb.g_L * memb.delta_T + I0 - Dist_VT
    
    if (Dist_VT<=0 or tau_ratio>=1 or memb.C<=0 or memb.g_L<=0 or memb.tau_w<=0 or memb.delta_T<=0):
        firing_rate=0
    else:
        dist_V_r = tau_ratio*(I0+memb.g_L*memb.delta_T*np.exp((V_r-memb.V_T)/memb.delta_T)-memb.g_L*(V_r-memb.E_L))  
        wV_V_r = -memb.g_L*(V_r-memb.E_L)+memb.g_L*memb.delta_T*np.exp((V_r-memb.V_T)/memb.delta_T)+I0
        w1= wV_V_r - dist_V_r
        w2= wV_V_r + dist_V_r

        delta_t1=0
        delta_t2=0
        if V_r >= memb.V_T:
            if w_r>=wV_V_r:
                w_ref=-memb.g_L*(memb.V_T-memb.E_L)+memb.g_L*memb.delta_T+I0+Dist_VT
                V_nullbound_inters=get_w_nullbound_inters(memb,I0,w_r,w_ref,1)            
                if len(V_nullbound_inters)<2:
                    delta_t1=quad(lambda x: ISI_integrate_part1and3(x,I0,memb,w_r),V_r,memb.V_T)[0]
                    delta_t2=0
                else:
                    V_cross_nullbound=np.min(V_nullbound_inters)
                    delta_t1=quad(lambda x: ISI_integrate_part1and3(x,I0,memb,w_r),V_r,V_cross_nullbound)[0]
                    delta_t2=quad(lambda x: ISI_integrate_part2(x,I0,memb),V_cross_nullbound,memb.V_T)[0]
                
                w_stop=w_end
                V1b=memb.V_T
                
            else:
                V1b=V_r
                w_stop=w_r
                   
        else:
            if(w_r < w2 and w_r > w1):
                delta_t2=quad(lambda x: ISI_integrate_part2(x,I0,memb),V_r,memb.V_T)[0]
                w_stop=w_end
            else:
                
                if (w_r <= w1):
                    signal=-1
                    w_ref=-memb.g_L*(memb.V_T-memb.E_L)+memb.g_L*memb.delta_T+I0-Dist_VT
                else:
                    signal=1
                    w_ref=-memb.g_L*(memb.V_T-memb.E_L)+memb.g_L*memb.delta_T+I0+Dist_VT # w_end
                
                V_nullbound_inters=get_w_nullbound_inters(memb,I0,w_r,w_ref,signal)
                
                if (not V_nullbound_inters or len(V_nullbound_inters)==1):
                    delta_t1=quad(lambda x: ISI_integrate_part1and3(x,I0,memb, w_r),V_r,memb.V_T)[0]
                    w_stop=w_r
                else:
                    V_bound=np.min(V_nullbound_inters)
                    delta_t1= quad(lambda x: ISI_integrate_part1and3(x,I0,memb, w_r),V_r,V_bound)[0]
                    delta_t2=quad(lambda x: ISI_integrate_part2(x,I0,memb),V_bound,memb.V_T)[0]
                    w_stop=w_end
                           
            V1b=memb.V_T

        if V1b>=memb.V_up:
            delta_t3=0
        else:
            
            delta_t3= quad(lambda x: ISI_integrate_part1and3(x,I0,memb,w_stop),V1b,memb.V_up)[0]

        ISI=np.asarray(delta_t1+delta_t2+delta_t3)
        firing_rate=1000/ISI

        return firing_rate


def Define_I_ref(I, memb):
       
    re = get_transFR(memb, I, 0, memb.V_r)

    if re is None:
        re = 0
        
    q = (re-200)**2
    
    return q

def latency_AdEx(memb, I):
    
    return quad(lambda V: memb.C/(I - memb.g_L*(V-memb.E_L) + memb.g_L*memb.delta_T*np.exp((V-memb.V_T)/memb.delta_T) ), memb.E_L, memb.V_T/2)[0]
    
def latency_LIF(memb, I):
    
    return memb.C*np.log(I/(I+memb.g_L*(memb.E_L-memb.V_T/2)))/memb.g_L
    

def random_param(N, par_mean, par_std, par_min, par_max, distr_flag):

    if par_std==0:
        if par_max==0:
            par = par_mean*np.ones(N)    
        else:
            par = par_min+(par_max-par_min)*np.random.random(size=N)       
    else:  
                         
        if distr_flag == 'normal':
            par = par_mean + par_std* np.random.normal(size=N)
            exc_ind = np.where((par<par_min) | (par>par_max))[0]

            par[exc_ind] =  par_min+(par_max-par_min)*np.random.random(size=len(exc_ind))
        elif distr_flag == 'uniform':
            par = par_min+(par_max-par_min)*np.random.random(size=N)
        elif distr_flag == 'log_normal':
            par = np.exp(np.random.normal(size=N) * par_std + par_mean)
            exc_ind = np.where((par<par_min) | (par>par_max))[0]
            par[exc_ind] = par_min+(par_max-par_min)*np.random.random(size=len(exc_ind))    
        
    return par


def get_commonneigh_recur(TS, Ncell):
    target, source = TS
    conn_mat = np.zeros((Ncell, Ncell))
    conn_mat[target, source] = 1
    conn_mat = conn_mat + conn_mat.transpose()
    conn_mat = conn_mat.astype(bool).astype(int)
    
    comneigh_mat = np.zeros((Ncell, Ncell))
    
    for tgt in range(Ncell):
        for src in range(tgt):
            comneigh_mat[tgt,src] = sum(conn_mat[tgt]*conn_mat[src]) - conn_mat[tgt, src]*conn_mat[tgt, tgt] - conn_mat[src, src]*conn_mat[src, tgt]

    comneigh_mat = comneigh_mat + comneigh_mat.transpose()
    
    return comneigh_mat.astype(int)    
    

def commonneighbour_report(TS_pairs, Ncells):
    Neigh_mat = get_commonneigh_recur(TS_pairs, Ncells)   
    N_neigh_connected = Neigh_mat[TS_pairs[0,:], TS_pairs[1, :]]  
    
    common_neigh = {}
    common_neigh['N'] = np.unique(Neigh_mat)
    common_neigh['all'] = np.zeros(len(common_neigh['N']))
    common_neigh['connected'] = np.zeros(len(common_neigh['N']))
    for i in range(len(common_neigh['N'])):
        common_neigh['all'][i] = np.sum((Neigh_mat == common_neigh['N'][i]).astype(int))
        common_neigh['connected'][i] = np.sum((N_neigh_connected==common_neigh['N'][i]).astype(int))
    
    common_neigh['ratio'] = common_neigh['connected']/common_neigh['all']

    return common_neigh
    

@dataclass
class Network(BaseClass):
    membr_params: xr.DataArray
    refractory_current: xr.DataArray
    syn_params: dict
    syn_pairs: np.ndarray
    group_distr: list
    seed: int or None
    basics: BaseClass
    
    def rheobase(self, neuron_idcs=None):
        g_L, V_T, E_L, delta_T = self.membr_params.loc[dict(param=['g_L', 'V_T', 'E_L', 'delta_T'])]
        Irheo =  self.basics.equations.rheobase(g_L, V_T, E_L, delta_T)
        
        if neuron_idcs is not None:
            Irheo = Irheo.loc[dict(cell_index=neuron_idcs)]
        
        return Irheo.values

    def save(self, path):
        if not os.path.isdir(path):
            os.mkdir(path)
        self.membr_params.to_netcdf(path+'//membr_params.nc')
        self.refractory_current.to_netcdf(path+'//refractory_current.nc')
        np.save(path+'//syn_pairs.npy', self.syn_pairs)
        os.mkdir(path+'//syn_params')
        for k in self.syn_params.keys():
            self.syn_params[k].to_netcdf(path+'//syn_params'+'//{}.nc'.format(k))
        
       
        stripe_list = []
        for stripe in self.group_distr:
            group_list = []
            for group in stripe:
                group_list.append(','.join(np.array(group).astype(str)))
            stripe_list.append(';'.join(group_list))
        group_distr = '\n'.join(stripe_list)      
        
        with open(path+'//group_distr.txt', 'w') as f:
            f.write(group_distr)
            
        with open(path+"//input.txt", 'w') as f:
            print(self.basics.struct.Ncells_prompt, self.basics.struct.stripes.N, sep='\n', end='', file=f)
        shutil.copyfile('BasicsSetup.py', path+'//BasicsSetup.py')
        
        if self.basics.scales is not None:
           with open(path+'//basics_scales.json', 'w') as f:
               json.dump(self.basics.scales, f)
        if self.seed is not None:
            with open(path+'//seed.txt', 'w') as f:
                f.write(str(self.seed))
                
           
    def load(path):
        
        print('REPORT: Loading Network from {}'.format(path), end='\n\n')
        
        refractory_current = xr.open_dataarray(path+'//refractory_current.nc')
        membr_params = xr.open_dataarray(path+'//membr_params.nc')
        channel = xr.open_dataarray(path+'//syn_params//channel.nc')
        delay = xr.open_dataarray(path+'//syn_params//delay.nc')
        spiking = xr.open_dataarray(path+'//syn_params//spiking.nc')
        STSP_params = xr.open_dataarray(path+'//syn_params//STSP_params.nc')
        syn_pairs = np.load(path+'//syn_pairs.npy')
            
        with open(path+'\\group_distr.txt', 'r') as f:
            group_distr_str = f.read()
        
        group_distr = []
        
        stripe_list = group_distr_str.split('\n')
        for stripe in stripe_list:
            group_distr.append([])
            group_list = stripe.split(';')
            for group in group_list:
                group_distr[-1].append([])
                cell_list = group.split(',')
                for cell in cell_list:
                    if len(cell_list)>1 or '' not in cell_list:
                        group_distr[-1][-1].append(int(cell))
        
        syn_params={'channel': channel, 'spiking': spiking,
                    'STSP_params': STSP_params, 'delay': delay}
        
        
        BasicsModule = import_module('{}.BasicsSetup'.format(path))
        
        with open(path+'\\input.txt', 'r') as f:
            prompt_str = f.read()
        Ncells_prompted, Nstripes = prompt_str.split('\n')
        Ncells_prompted = int(Ncells_prompted) 
        Nstripes = int(Nstripes)
        
        basics_scales = None
        if 'basics_scales.json' in os.listdir(path):
            with open(path+'\\basics_scales.json', 'r') as f:
                basics_scales = json.load(f)
    
        seed = None
        if 'seed.txt'in os.listdir(path):
            with open(path+'//seed.txt', 'r') as f:
                seed = int(f.read())
                
        basics = BasicsModule.basics_setup(Ncells_prompted, Nstripes, basics_scales, disp=False)
        
        
        
        print('REPORT: Network loaded\n')
        
        return Network(membr_params, refractory_current, syn_params, syn_pairs, group_distr, seed, basics)
         
    def group_idcs(self, group):
        if isinstance(group, str):
            names = []
            for gr in self.basics.struct.groups.sets[group]:
                names.append(self.basics.struct.groups.idcs[gr])
            return names
        elif isinstance(group, int):
            return [group]
        
    def neuron_idcs(self, groupstripe_list):
        
        if isinstance(groupstripe_list[0], (str,int)):
            groupstripe_list = [groupstripe_list]
            
        neuron_idcs = []
        for groupstripe in groupstripe_list:
            group, stripe_idc = groupstripe
            group_idc = self.group_idcs(group)
            for idc in group_idc:
                neuron_idcs.extend(self.group_distr[stripe_idc][idc])
        
        return np.array(neuron_idcs)
    
    def syn_idcs_from_neurons(self, target, source, channel=None):
        
        syn_idcs = np.arange(self.syn_pairs.shape[1])
        target_isin = np.isin(self.syn_pairs[0,:], np.array(target))
        source_isin = np.isin(self.syn_pairs[1,:], np.array(source))
        both_isin = target_isin & source_isin
        
        syn_idcs = syn_idcs[both_isin]
        
        if channel is not None:
            ch_isin = []
            if isinstance(channel, str):
                channel = [channel]
            for ch in channel:
                ch_idc = np.where(self.syn_params['channel'].values[0] == self.basics.syn.channels.names.index(ch))
                ch_isin.append(syn_idcs[np.isin(syn_idcs, ch_idc)])
            ch_isin = np.concatenate(ch_isin) 
            syn_idcs = syn_idcs[np.isin(syn_idcs, ch_isin)]
            
        return syn_idcs
        
    
    def syn_idcs_from_groups(self, target_groupstripe_list, source_groupstripe_list, channel=None):
        
        if isinstance(target_groupstripe_list[0], (str,int)):
            target_groupstripe_list = [target_groupstripe_list]
        
        if isinstance(source_groupstripe_list[0], str) or isinstance(source_groupstripe_list[0], int):
            source_groupstripe_list = [source_groupstripe_list]
        
        syn_idcs = np.array([])
        all_syn = np.arange(self.syn_pairs.shape[1])
        for target_groupstripe in target_groupstripe_list:
            for source_groupstripe in source_groupstripe_list:
                targetgroup, targetstripe_idc = target_groupstripe
                sourcegroup, sourcestripe_idc = source_groupstripe
                
                targetgroup_idc = self.group_idcs(targetgroup)
                sourcegroup_idc = self.group_idcs(sourcegroup)
                
                targetneurons_idc = []
                for gr in targetgroup_idc:
                    targetneurons_idc.extend(self.group_distr[targetstripe_idc][gr])
               
                sourceneurons_idc = []
                for gr in sourcegroup_idc:
                    sourceneurons_idc.extend(self.group_distr[sourcestripe_idc][gr])
                
                target_isin = np.isin(self.syn_pairs[0,:], np.array(targetneurons_idc))
                source_isin = np.isin(self.syn_pairs[1,:], np.array(sourceneurons_idc))
                both_isin = target_isin & source_isin
                syn_idcs = np.concatenate((syn_idcs, all_syn[both_isin]))
        
        syn_idcs = syn_idcs.astype(int)   
        
        if channel is not None:
            ch_isin = []
            if isinstance(channel, str):
                channel = [channel]
            for ch in channel:
                ch_idc = np.where(self.syn_params['channel'].values[0] == self.basics.syn.channels.names.index(ch))
                ch_isin.append(syn_idcs[np.isin(syn_idcs, ch_idc)])
            ch_isin = np.concatenate(ch_isin) 
            syn_idcs = syn_idcs[np.isin(syn_idcs, ch_isin)]
        return syn_idcs
    
    
    
    
    
if __name__ == '__main__':
    
    basics_scales={}
    # gmax_scales = [(dict(target=GROUPS_SETS['PC'], source=GROUPS_SETS['PC']), 2)]    
    # basics_scales['gmax_mean'] = gmax_scales
    network = network_setup(1000, 1)
    # network.save('New6')
    # network = Network.load('New5')

    
    
    