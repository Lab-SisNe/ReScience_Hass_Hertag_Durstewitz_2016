import numpy as np
from ParamSetup import ParamSetup
from scipy.integrate import quad
from scipy.optimize import fmin
from scipy.optimize import fsolve 
from itertools import product, combinations


def NetworkSetup(N1, Nstripes, scales, clustering):   

    NumSynTypes = 2
    
    (NTypes, CellPar, V0Par, k_trans, ParCov, N_sig, N_max, N_min,
    SynTypes, Syngmax, Syndelay, Synpfail, STSP_prob, STSP_types, pCon, cluster_flag,
    S_sig, S_max, S_min, p_fail,
    STSP_E1, STSP_E2, STSP_E3, STSP_I1, STSP_I2, STSP_I3,
    STSP_setup, get_SynType, get_STSP_prob, get_STSP_types, STypPar,
    ConParStripes,
    ) = ParamSetup(N1, Nstripes, scales)
  
    N = int(round(np.sum(NTypes), 0))             
    NTypesN = len(NTypes)       

    NeuPar = np.zeros((len(CellPar[:,0])+2,Nstripes*N))
    V0 = np.zeros((2,Nstripes*N))
    SPMtx=np.zeros((Nstripes*(N),Nstripes*(N), NumSynTypes))
    
    NN = 0
    SynPar=[]
    Nsyn=0
    delay_back = np.size([])
    
    ind_NeuPar = [i for i in range(0, 9)]
    valid_ind_NeuPar = [0, 1, 2, 3, 4, 5, 7, 8, 9]
    
    recur_clustering = clustering
    
    # -------------- Set neuron parameters and most synaptic connections -----------------
 
    group_distr = [ [[] for j in range(Nstripes)] for i in range (NTypesN)]

    t_lat = np.zeros(N*Nstripes)
    t_lat_LIF = np.zeros(N*Nstripes)
    adapt = np.zeros(N*Nstripes)  
  
    for ii in range(Nstripes):   
        for i in range(NTypesN):
            iset = (np.asarray([k for k in range(int(round(NTypes[i], 0)))])+NN).astype(int)
            group_distr[i][ii].extend(iset)
            ind_out = iset
        
            iii=0
    
            while ind_out.size and iii < 999:
                NeuPar_multi=np.random.multivariate_normal(np.zeros(len(ParCov[i][:,0])),ParCov[i],len(ind_out)).transpose()

                for j in ind_NeuPar:
                    _ = inv_transform_distribution2(NeuPar_multi[j,:],k_trans[j,i],CellPar[j+1,i],N_sig[j+1,i],N_min[valid_ind_NeuPar[j]+1, i])                  
                    NeuPar_multi[j,:] = _
                    
                for j in ind_NeuPar:
                    NeuPar[valid_ind_NeuPar[j],ind_out] = NeuPar_multi[j,:]
                                   
                NeuPar[0,ind_out] = NeuPar[0,ind_out]*NeuPar[1,ind_out]   # C from gL and tau
                ind_out_para = np.asarray([])
                
                for j in ind_NeuPar:
                    bound_idc = np.where((NeuPar[valid_ind_NeuPar[j],iset]<N_min[valid_ind_NeuPar[j]+1,i]) | (NeuPar[valid_ind_NeuPar[j],iset] > N_max[valid_ind_NeuPar[j]+1,i]))[0]
                    ind_out_para = np.concatenate((ind_out_para, bound_idc))
                    
                ind_out_V2 = np.where((NeuPar[8,iset]>=NeuPar[9,iset]))[0]   # Vr must be smaller than Vth
                ind_out_tau = np.where((NeuPar[0,iset]/NeuPar[1,iset]<N_min[11,i]) |( NeuPar[0,iset]/NeuPar[1,iset]>N_max[11,i]))[0] # check tau boundaries)
                ind_out_tcw = np.where((NeuPar[5,iset]<=NeuPar[0,iset]/NeuPar[1,iset]))[0]
                ind_out_nan = np.argwhere(np.isnan(NeuPar[:, iset]))[:, 1]
                all_idc = np.concatenate((ind_out_para,ind_out_V2,ind_out_tau , ind_out_tcw, ind_out_nan ))

                ind_out=iset[np.unique(all_idc).astype(int)] 

                iii=iii+1
                
            iii=0

            while ind_out.size and iii<1000:
                if iii==0:
                    print('WARNING: Number of trials for the full multivariate distribution has been exceeded. Use multivariate uniform distribution for the rest.\n');
                
                NeuPar_multi=np.random.uniform(N_min[np.asarray(valid_ind_NeuPar)+1,i], N_max[np.asarray(valid_ind_NeuPar)+1,i], len(ind_out)).transpose()
                
                for j in ind_NeuPar:
                    NeuPar[valid_ind_NeuPar[j],ind_out] = NeuPar_multi[j,:]
                
                NeuPar[0,ind_out] = NeuPar[0,ind_out]*NeuPar[1,ind_out]                # C from gL and tau
                NeuPar[6,ind_out] = np.zeros((np.size(NeuPar[3,ind_out])))                    # a=0
                
                ind_out_para = np.asarray([])
                
                for j in ind_NeuPar:
                    bound_idc = np.where((NeuPar[valid_ind_NeuPar[j],iset]<N_min[valid_ind_NeuPar[j]+1,i]) | (NeuPar[valid_ind_NeuPar[j],iset] > N_max[valid_ind_NeuPar[j]+1,i]))[0]
                    ind_out_para = np.concatenate((ind_out_para, bound_idc))
                    
                ind_out_V2 = np.where((NeuPar[8,iset]>=NeuPar[9,iset]))[0]   # Vr must be smaller than Vth
                ind_out_tau = np.where((NeuPar[0,iset]/NeuPar[1,iset]<N_min[11,i]) |( NeuPar[0,iset]/NeuPar[1,iset]>N_max[11,i]))[0] # check tau boundaries)
                ind_out_tcw = np.where((NeuPar[5,iset]<=NeuPar[0,iset]/NeuPar[1,iset]))[0]
                         
            if ind_out.size:
                print('ERROR: Number of trials for univariate uniform distribution has been exceeded.\n')
            
            t_lat_act = np.zeros(np.size(iset));
            t_lat_LIF_act = np.zeros(np.size(iset));
            adapt_act = np.zeros(np.size(iset));
            
            for j in range(len(iset)):

                I   = 500
                C = NeuPar[0,iset[j]]
                gL = NeuPar[1,iset[j]]
                EL = NeuPar[2,iset[j]] 
                sf = NeuPar[3,iset[j]] 
                Vup = NeuPar[4, iset[j]]
                tcw = NeuPar[5,iset[j]]                 
                b = NeuPar[7,iset[j]]
                Vr = NeuPar[8,iset[j]]
                Vth = NeuPar[9,iset[j]]
                Irheo = gL*(Vth-EL) -gL * sf
                Iref = fmin(Define_I_ref, Irheo +100, disp=False, args=(NeuPar[:,iset[j]],))
                NeuPar[10, iset[j]] = Iref
                NeuPar[11, iset[j]] = Vr                                                      

                t_lat_act[j] = quad(lambda V: C/(I - gL*(V-EL) + gL*sf*np.exp((V-Vth)/sf) ), EL, Vth/2)[0]
                t_lat_LIF_act[j] = C*np.log(I/(I+gL*(EL-Vth/2)))/gL

                I = [i for i in range(25,301,25)]
                f1 = np.zeros(len(I))
                f  = np.zeros(len(I))
                
                for k in range(len(I)):
                    f1[k] = FRsimpAdEx(NeuPar[:,iset[j]],I[k],0,NeuPar[2,iset[j]],[])
                    f[k]  = FRsimpAdEx(NeuPar[:,iset[j]],I[k],[],[],[])
                
                idc0 = np.where((f>0))[0]
                
                if np.isnan(np.median(f1[idc0] / f[idc0])):
                    adapt_act[j] = 0
                else:
                    adapt_act[j] = np.median(f1[idc0] / f[idc0])
          
            t_lat[iset] = t_lat_act
            t_lat_LIF[iset] = t_lat_LIF_act
            adapt[iset] = adapt_act            
   
            V0[1,iset]=np.zeros((1,len(iset)))
            
            for j in range(len(iset)):
                V0[0,iset[j]]=NeuPar[2,iset[j]]  
            
            NN = NN+NTypes[i]
        
        # Redistribute neuron types: I-L -> I-L-d          
        
        ind_L23 = np.asarray(group_distr[1][ii] +  group_distr[2][ii])
        
        if len(ind_L23 > 0):
            idc_0 = np.where((t_lat[ind_L23] - t_lat_LIF[ind_L23]>0))[0]    
            ind_L23_Ld = ind_L23[idc_0]
            ind_L23_L = np.setdiff1d(ind_L23, ind_L23_Ld)
            NTypes[1] = len(ind_L23_L)
            NTypes[2] = len(ind_L23_Ld)
            group_distr[1][ii] = list(ind_L23_L)
            group_distr[2][ii] = list(ind_L23_Ld)
        
        ind_L5  = np.asarray(group_distr[8][ii] +  group_distr[9][ii])
        
        if len(ind_L5)>0:
            idc_0 = np.where((t_lat[ind_L5] - t_lat_LIF[ind_L5]>0))[0]           
            ind_L5_Ld = ind_L5[idc_0]
            ind_L5_L = np.setdiff1d(ind_L5, ind_L5_Ld)           
            NTypes[8] = len(ind_L5_L)
            NTypes[9] = len(ind_L5_Ld)
            group_distr[8][ii] = list(ind_L5_L)
            group_distr[9][ii] = list(ind_L5_Ld)
            
        # Redistribute neuron types: I-CL -> I-CL-AC
        
        ind_L23 = np.asarray(group_distr[3][ii] + group_distr[4][ii])
        
        if len(ind_L23)>0:
            idc_0 = np.where((adapt[ind_L23]>1.5834))[0]
            ind_L23_Ld = ind_L23[idc_0]
            ind_L23_L = np.setdiff1d(ind_L23, ind_L23_Ld)    
            NTypes[3] = len(ind_L23_L)
            NTypes[4] = len(ind_L23_Ld)
            group_distr[3][ii] = list(ind_L23_L)
            group_distr[4][ii] = list(ind_L23_Ld)
        
        ind_L5  = np.asarray(group_distr[10][ii] + group_distr[11][ii])
        
        if len(ind_L5)>0:
            idc_0 = np.where((adapt[ind_L5]>1.5834))[0]
            ind_L5_Ld = ind_L5[idc_0]
            ind_L5_L = np.setdiff1d(ind_L5, ind_L5_Ld)
            NTypes[10] = len(ind_L5_L)
            NTypes[11] = len(ind_L5_Ld)
            group_distr[10][ii] = list(ind_L5_L)
            group_distr[11][ii] = list(ind_L5_Ld)

        # intra-celltype connections
        
        for i in range(NTypesN):
            if len(group_distr[i][ii]):    
                if cluster_flag[i,i] == 1 and recur_clustering:
                    X, idc = SetCon_CommonNeighbour_Recur(Nsyn, len(group_distr[i][ii]) , pCon[i, i], 0.47)
                else:
                    X,idc = SetCon(Nsyn, len(group_distr[i][ii]), len(group_distr[i][ii]), pCon[i, i])         # ... or without common neighbour rule

                if len(idc)>0:
                    SynTypes_act = get_SynType[SynTypes[i, i]]
                
                    for k in range(len(SynTypes_act)):    
                        ST = SynTypes_act[k]
                        _ = X + (k) * len(idc)
                        cart_prod = list(product(group_distr[i][ ii], group_distr[i][ ii]))
                        Idc1, Idc2 = map(list, zip(*cart_prod))
                        SPMtx[Idc1, Idc2, k] = _.flatten()              
                        idc_back = idc - 1
                        j = i
                        
                        SynPar, Nsyn, idc = SetSyn(SynPar, Nsyn, idc, ST, Syngmax[i, i], Syndelay[i, i], Synpfail[i, i], STSP_types[i, i], STSP_prob[i, i], get_STSP_prob, get_STSP_types,STSP_setup,
                                                    S_sig[i, i, :], S_max[i, i, :], S_min[i, i,:])
                        
                        if k ==0:
                            delay_back = SynPar[5, idc_back]
                        else:
                            SynPar[5, idc_back] = delay_back
                            delay_back = []
            
        for i in range(NTypesN):                              # loop over output neurons
            if len(group_distr[i][ii]):    
                for j in np.setdiff1d(np.arange(NTypesN),[i,]) :             # loop over input neurons (except output neuron)
                    if len(group_distr[j][ii]):                    
                                          
                        X, idc = SetCon(Nsyn, len(group_distr[i][ii]), len(group_distr[j][ii]), pCon[i,j])                                                                             
                                               
                        if np.size(idc):
                            SynTypes_act = SynTypes_act = get_SynType[SynTypes[i, j]]
                        
                            for k in range(len(SynTypes_act)):                                                                                                                          
                                ST = SynTypes_act[k]   
                                _ = X + (k) * len(idc)
                                cart_prod = list(product(group_distr[i][ ii], group_distr[j][ ii]))
                                Idc1, Idc2 = map(list, zip(*cart_prod))
                                SPMtx[Idc1, Idc2, k] = _.flatten()           
                                idc_back = idc - 1

                                SynPar,Nsyn,idc = SetSyn(SynPar,Nsyn,idc,ST, Syngmax[i, j], Syndelay[i, j], Synpfail[i, j], STSP_types[i, j], STSP_prob[i, j],get_STSP_prob, get_STSP_types,STSP_setup,
                                                           S_sig[i,j,:], S_max[i,j,:2], S_min[i,j,:]) # set synapse parameters / +2 replaced by +1
 
                                if k==0:
                                    delay_back = SynPar[5,idc_back]
                                else:
                                    SynPar[5,idc_back] = delay_back
                                    delay_back = []
                       
    print('REPORT: Synaptic connectivity and parameters generated\n')
    
    # # ------------------------- Define inter-stripe connections ---------------------------------
    
    if Nstripes>1:
        pConStripes, coefStripes, interStripes = ConParStripes
        
        for i in range(len(pConStripes)):
            i_act = pConStripes[i][0]
            j_act = pConStripes[i][1]
            
            for j in range(Nstripes):             
                for kk in range(len(interStripes[i])):
                    target_act = j + interStripes[i][kk]
                    
                    while target_act<0:
                        target_act = target_act + Nstripes
                    
                    while target_act>=Nstripes:
                        target_act = target_act - Nstripes
                    
                    dist_act = abs(interStripes[i][kk])
                    p_act = pCon[i_act,j_act]*np.exp(-dist_act/coefStripes[i][0])

                   

                    
                    X,idc = SetCon(Nsyn, len(group_distr[i_act][target_act]), len(group_distr[j_act][j]), p_act)                                                                                                     # ... or without common neighbour rule
                    
                    SynTypes_act = get_SynType[SynTypes[i_act, j_act]]
                    
                    for k in range(len(SynTypes_act)):                                                   
                                                                   
                        ST = SynTypes_act[k]                                                                  
                        _ = (X)+(k)*len(idc)
                        
                        cart_prod = list(product(group_distr[i_act][target_act], group_distr[j_act][j]))
                        Idc1, Idc2 = map(list, zip(*cart_prod))
                        SPMtx[Idc1, Idc2, k] = _.flatten()

                        gmax_act = Syngmax[i_act, j_act]*np.exp(-dist_act/coefStripes[i][0]) 
                        delay_act = Syndelay[i_act, j_act] + coefStripes[i][1]*dist_act
                        S_sig_act = [0, 0]
                        S_sig_act[0] = S_sig[i_act,j_act,0]*np.exp(-dist_act/coefStripes[i][0])
                        S_sig_act[1] = S_sig[i_act,j_act,1] + coefStripes[i][1]*dist_act
                        
                        idc_back = idc - 1
    
                        p_fail_act = 0.3
                        
                        SynPar,Nsyn,idc =SetSyn(SynPar,Nsyn,idc,ST, gmax_act, delay_act, p_fail_act,
                                                STSP_types[i_act, j_act], STSP_prob[i_act, j_act], get_STSP_prob, get_STSP_types,STSP_setup,
                                                S_sig_act, S_max[i_act,j_act,:], S_min[i_act,j_act,:])
    
                        if k==0:
                            delay_back = SynPar[5,idc_back]
                        else:
                            SynPar[5,idc_back] = delay_back
                            delay_back = []
    
    return NeuPar, V0, STypPar, SynPar,SPMtx, group_distr

    
def SetCon(Nsyn, Nin, Nout, pCon):
        
    k1 = np.arange(Nin*Nout)
    k = k1[:int(round(pCon*Nin*Nout, 0))]
    np.random.shuffle(k1)
    X = np.zeros((Nin, Nout))
    idc = Nsyn + np.arange(len(k)) + 1
   
    X_ = X.transpose().flatten() 
    X_[k] = idc
    X = np.reshape(X_, X.transpose().shape).transpose()

    return X, idc


def SetSyn(SynPar,Nsyn,idc,ST,gmax, delay, p_fail, STSP_types, STSP_prob, get_STSP_prob, get_STSP_types, STSP_setup, S_sig,S_max,S_min):
       
    idc_SynPar = idc - 1
    
    if not np.size(SynPar):
        SynPar = np.zeros((7, max(idc_SynPar)+1))

    if np.max(idc)> SynPar.shape[1]:
        SynPar_temp = SynPar
        SynPar = np.zeros((7, max(idc_SynPar)+1))
        SynPar[:, :SynPar_temp.shape[1]] = SynPar_temp
    
    SynPar[0, idc_SynPar] = ST
    
    mean_gmax = np.log((gmax**2)/np.sqrt(S_sig[0]**2 + gmax **2))
    std_gmax = np.sqrt(np.log(S_sig[0]**2/gmax**2 + 1))
    
    if gmax == 0 or S_sig[0] ==0:
        SynPar[4, idc_SynPar] = rand_par(len(idc_SynPar), gmax, S_sig[0], S_min[0]*gmax, S_max[0] * gmax, 0)
    else:
        SynPar[4, idc_SynPar] = rand_par(len(idc_SynPar), mean_gmax, std_gmax, S_min[0]*gmax, S_max[0]*gmax, 2)
    
    SynPar[5, idc_SynPar] = rand_par(len(idc_SynPar), delay, S_sig[1], S_min[1]*delay, S_max[1]*delay, 0)
    SynPar[6, idc_SynPar] = p_fail

    ConParSTSP0_arr = np.asarray(get_STSP_prob[STSP_prob])

    if np.size(idc_SynPar):    
        idc_frac = (np.round(ConParSTSP0_arr * len(idc_SynPar))).astype(int)

        if sum(idc_frac) > len(idc_SynPar):
            ind = np.argmin(ConParSTSP0_arr*len(idc_SynPar) - np.floor(ConParSTSP0_arr*len(idc_SynPar)))
            idc_frac[ind] = idc_frac[ind] - 1           
        elif sum(idc_frac) < len(idc_SynPar):
                ind = np.argmax(idc_frac - np.floor(idc_frac))
                idc_frac[ind] = idc_frac[ind] + 1

    idc_act = (Nsyn + np.arange(idc_frac[0])).astype(int)
    Nsyn_act = Nsyn + idc_frac[0]
    
    STSP = get_STSP_types[STSP_types]
    
    for i in range(len(ConParSTSP0_arr)):
        STSP_params = STSP_setup[STSP[i]]
        use = STSP_params[0, 0]
        tc_rec = STSP_params[0, 1]
        tc_fac = STSP_params[0, 2]
        std_use    = STSP_params[1, 0]
        std_tc_rec = STSP_params[1, 1]
        std_tc_fac = STSP_params[1, 2]
        
        SynPar[1, idc_act] =    rand_par(len(idc_act), use,    std_use,    0,    1, 0)     
        SynPar[2,idc_act] = rand_par(len(idc_act), tc_rec, std_tc_rec, 0, 1500, 0)
        SynPar[3,idc_act] = rand_par(len(idc_act), tc_fac, std_tc_fac, 0, 1500, 0)
    
        if i < len(ConParSTSP0_arr)-1 and len(ConParSTSP0_arr)>1:
            idc_act = (Nsyn_act + np.arange(idc_frac[i+1])).astype(int) # +1 taken out
            Nsyn_act = Nsyn_act + idc_frac[i+1]
    
    Nsyn += len(idc)
    idc = Nsyn+np.arange(len(idc)  ) + 1
        
    return SynPar,Nsyn,idc

def sourcetarget(SPMtx):
    
    target_arr = np.zeros((np.max(SPMtx).astype(int)))
    source_arr = np.zeros((np.max(SPMtx).astype(int)))
    
    for i in range(SPMtx.shape[0]): # target
        for j in range(SPMtx.shape[1]): # source
            _ = SPMtx[i, j, 0].astype(int)
            if _ > 0:
                _ -= 1
                target_arr[_] = i
                source_arr[_] = j
                _ = SPMtx[i, j, 1].astype(int)
                if _ > 0:
                    _ -= 1
                    target_arr[_] = i
                    source_arr[_] = j
                    
    return target_arr.astype(int), source_arr.astype(int)
                            
def inv_transform_distribution2(X_trans, k, mean_X, std_X, min_X):
    
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


def FRsimpAdEx(values,I,w0=[],V0=[],bmin=[]):

    Cm,gL,EL,sf,Vup,tcw,a,b,Vr, Vth=values[:10]

    tau=Cm/gL
    f=tau/tcw
    X_Vth=f*(I+gL*sf-gL*(Vth-EL))                  
    w_end=-gL*(Vth-EL)+gL*sf+I-X_Vth

    if not np.size(w0):
        if b!=0:
            w_r=w_end+b
        else:
            w_r=0    
    else:
        w_r=w0
    
    if not np.size(bmin) or bmin<0:
        bmin=0;
    
    if not np.size(V0):
        V_r=Vr
    else:
        V_r=V0

    if (X_Vth<=0 or f>=1 or b<bmin or Cm<=0 or gL<=0 or tcw<=0 or sf<=0):
        fr=0
    else:
        X_Vr=f*(I+gL*sf*np.exp((V_r-Vth)/sf)-gL*(V_r-EL))  
        wV_Vr=-gL*(V_r-EL)+gL*sf*np.exp((V_r-Vth)/sf)+I
        w1=wV_Vr-X_Vr
        w2=wV_Vr+X_Vr

        t1=0
        t2=0
        if V_r >= Vth:
            if w_r>=wV_Vr:
                w_ref=-gL*(Vth-EL)+gL*sf+I+X_Vth
                Vfix=IntPoints_simpAdEx(values,I,w_r,w_ref,1)
                i=0
                j=w_r
                k=0               
                if (not Vfix or len(Vfix)==1):
                    t1=quad(lambda x: ISI_int_piecewise(x,I,values,i,j,k),V_r,Vth)[0]
                    t2=0
                else:
                    Vlb=np.min(Vfix)
                    t1=quad(lambda x: ISI_int_piecewise(x,I,values,i,j,k),V_r,Vlb)[0]
                    m=f*gL-gL
                    i=m
                    j=(1-f)*I-m*EL
                    k=(1-f)*gL*sf
                    t2=quad(lambda x: ISI_int_piecewise(x,I,values,i,j,k),Vlb,Vth)[0]
                
                w_stop=w_end
                Vlb=Vth
                
            else:
                Vlb=V_r
                w_stop=w_r
                   
        else:
            if(w_r < w2 and w_r > w1):
                m=f*gL-gL;
                i=m
                j=(1-f)*I-m*EL
                k=(1-f)*gL*sf
                t2=quad(lambda x: ISI_int_piecewise(x,I,values,i,j,k),V_r,Vth)[0]
                w_stop=w_end
            else:
                i=0
                j=w_r
                k=0
                
                if (w_r <= w1):
                    ns=-1
                    w_ref=-gL*(Vth-EL)+gL*sf+I-X_Vth
                else:
                    ns=1
                    w_ref=-gL*(Vth-EL)+gL*sf+I+X_Vth # w_end
                
                Vfix=IntPoints_simpAdEx(values,I,w_r,w_ref,ns)
                
                if (not Vfix or len(Vfix)==1):
                    t1=quad(lambda x: ISI_int_piecewise(x,I,values,i,j,k),V_r,Vth)[0]
                    w_stop=w_r
                else:
                    Vlb=np.min(Vfix)
                    t1= quad(lambda x: ISI_int_piecewise(x,I,values,i,j,k),V_r,Vlb)[0]
                    m=f*gL-gL
                    i=m
                    j=(1-f)*I-m*EL
                    k=(1-f)*gL*sf
                    t2=quad(lambda x: ISI_int_piecewise(x,I,values,i,j,k),Vlb,Vth)[0]
                    w_stop=w_end
                           
            Vlb=Vth

        if Vlb>=Vup:
            t3=0
        else:
            i=0
            j=w_stop
            k=0
            t3= quad(lambda x: ISI_int_piecewise(x,I,values,i,j,k),Vlb,Vup)[0]

        ISI=np.asarray(t1+t2+t3)
        fr=1000/ISI

        return fr


def ISI_int_piecewise(V,I,values,i,j,k):
    
    Cm,gL,EL,sf,Vup,tcw,a,b,Vr,Vth = values[:10]

    F=(1/Cm)*(I-(i*V+j+k*np.exp((V-Vth)/sf))+gL*sf*np.exp((V-Vth)/sf)-gL*(V-EL))
    f=1/F
    
    return f


def IntPoints_simpAdEx(values,I,w0,w_ref,i):

    Cm,gL,EL,sf,Vup,tcw,a,b,Vr,Vth= values[:10]
    numcor = 0
    
    Vfix=[None, None]
    
    if w0>w_ref:
        f=Cm/(gL*tcw)
        G= lambda V: ((1+i*f)*(I-gL*(V-EL)+gL*sf*np.exp((V-Vth)/sf)) - w0)

        lb=EL+(I-(w0/(1+i*f)))/gL-0.1
        ub=EL+sf+(I-(w0/(1+i*f)))/gL
        
        while (np.sign(G(ub))/np.sign(G(lb))==1):
            lb=EL+(I-(w0/(1+i*f)))/gL-numcor
            ub=EL+sf+(I-(w0/(1+i*f)))/gL
            numcor=numcor+1
            
            if(numcor>1000):
                print('Error in numcor')
                
        Vfix[0] =fsolve(G,(lb+ub)/2)
        lb=EL+sf+(I-(w0/(1+i*f)))/gL
        ub=Vup
        numcor=0
        
        while(np.sign(G(ub))/np.sign(G(lb))==1):
            lb=EL+sf+(I-(w0/(1+i*f)))/gL-numcor
            ub=Vup
            numcor=numcor+1
            if(numcor>1000):
                print('Error')
        
        Vfix[1]= fsolve(G,(lb+ ub)/2)      
    elif w0==w_ref:
        Vfix=Vth      
    else:
        Vfix=[] 
    
    return Vfix


def Define_I_ref(I0, Par):
    
    re = FRsimpAdEx(Par, I0, 0, [], [])

    if type(re) == type(None):
        re = 0
        
    q = (re-200)**2
    
    return q


def rand_par(N, par_mean, par_std, par_min, par_max, distr_flag):

    q = 0
    partial = 0
    exc = 0
    if par_std==0:
        if par_max==0:
            par = par_mean*np.ones(N)    
        else:
            par = par_min+(par_max-par_min)*np.random.random(size=N)       
    else:                        
        if distr_flag == 0:
            par = par_mean + par_std* np.random.normal(size=N)
            exc_ind = np.where((par<par_min) | (par>par_max))[0]

            par[exc_ind] =  par_min+(par_max-par_min)*np.random.random(size=len(exc_ind))
        elif distr_flag == 1:
            par = par_min+(par_max-par_min)*np.random.random(size=N)
        elif distr_flag == 2:
            par = np.exp(np.random.normal(size=N) * par_std + par_mean)
            exc_ind = np.where((par<par_min) | (par>par_max))[0]
            par[exc_ind] = par_min+(par_max-par_min)*np.random.random(size=len(exc_ind))            
        
    return par


def SetCon_CommonNeighbour_Recur(Nsyn, Npop, pCon, pSelfCon):
    
    slope = 20*3.9991 / Npop
    
    X, _ = SetCon(Nsyn, Npop, Npop, pCon)
    N_neigh = find_neigh_recur(X)
    p0 = p_calc_recur(X, N_neigh, pCon, pSelfCon, slope)
    
    p = np.zeros((p0.shape[0], 2))
    p[:, 0] = p0[:,0]
    p[:, 1] = p0[:,5]
    
    pair_id_selected = []
    
    X_nondiag = X.copy()
    
    for i in range(X.shape[0]):
        X_nondiag[i,i] = 0
        
    X_diag = X.copy()
    
    for i in range(X.shape[0]):
        X_diag[i,i] = 1
           
    k = set(map(tuple, zip(*np.where(X_nondiag>0))))
    
    X_int = (X_nondiag>0).astype(int)
    X_recur = X_int * X_int.transpose()
    X_recur_tril = np.tril(X_recur)
    X_none = (X_diag + X_diag.transpose()==0).astype(int) 

    k_recur_tril = set(map(tuple, zip(*np.where(X_recur_tril>0))))
    k_recur = set(map(tuple, zip(*np.where(X_recur>0))))
    k_missing = set(map(tuple, zip(*np.where(X_none>0))))
    
    pairs_select = []
    
    for i in range(p.shape[0]):
        pair_id_act = set(map(tuple, zip(*np.where(N_neigh==p[i, 0]))))     
        k_missing_act = k_missing.intersection(pair_id_act)
        pairs_id_old = k.intersection(pair_id_act)
        
        pairs_id_recur_old = k_recur_tril.intersection(pair_id_act)
        pairs_id_recur_old2 = k_recur.intersection(pair_id_act)
        pairs_id_uni_old = pairs_id_old.difference(pairs_id_recur_old2)
        
        N_rec = (np.floor(p[i,1]*pSelfCon/2 * len(pair_id_act))).astype(int)
        N_uni = (np.ceil(p[i,1] * len(pair_id_act)) - 2 * N_rec).astype(int)

        N_old_recur = np.sum(X_recur)
        N_old_uni = np.sum(X_int) - N_old_recur
        
        if len(pairs_id_recur_old) >= N_rec:
            k1_old = np.arange(len(pairs_id_recur_old))
            np.random.shuffle(k1_old)
            k1_select = np.asarray(k1_old[:N_rec])
            
            for cell1, cell2 in np.asarray(list(pairs_id_recur_old))[k1_select]:          
                pairs_select.extend([(cell1, cell2), (cell2, cell1)])
          
        else:
            for cell1, cell2 in list(pairs_id_recur_old):
                pairs_select.extend([(cell1, cell2), (cell2, cell1)])
           
            N_new = N_rec  - len(pairs_id_recur_old)
            missing_pairs_id = list(k_missing_act)
            k1_new = np.arange(len(missing_pairs_id))
            np.random.shuffle(k1_new)
            k1_select = np.asarray(k1_new[:N_new])
            select = []
            
            kk = 0
            kkk = 0
            
            while kk < N_new:
                cell1, cell2 = missing_pairs_id[k1_new[kkk]]
                kkk+=1
                
                if cell1>cell2:
                    pairs_select.extend([(cell1, cell2), (cell2, cell1)])
                    select.extend([(cell1, cell2), (cell2, cell1)])
                    kk += 1
       
            if len(pairs_select):
                _ = np.unique(np.vstack(pairs_select), axis=0, return_counts=1)[1]

            k_missing_act = k_missing_act.difference(set(select))
        
        if len(pairs_id_uni_old) > N_uni:            
            k1_old = np.arange(len(pairs_id_uni_old))
            np.random.shuffle(k1_old)
            k1_select = np.asarray(k1_old[:N_uni])
            pairs_select.extend(np.asarray(list(pairs_id_uni_old))[k1_select])        
        else:
            pairs_select.extend(list(pairs_id_uni_old))
            N_new = N_uni - len(pairs_id_uni_old)

            missing_pairs_id = np.asarray(list(k_missing_act))
            k1_new = np.arange((len(missing_pairs_id)))
            np.random.shuffle(k1_new)

            k1_select = np.asarray(k1_new[:N_new])
            pairs_select.extend(missing_pairs_id[k1_select])
    
    X_diag = X * np.identity(X.shape[0])

    pairs_diag = [(i, i) for i in range(X.shape[0])]
    pairs_diag_old = list(map(tuple, zip(*np.where(X_diag>0))))    
    N_diag_old = np.sum((X_diag>0).astype(int))
    
    if N_diag_old > int(round(pSelfCon*len(X_diag), 0)):
        k1_old = np.arange(N_diag_old)
        np.random.shuffle(k1_old)
        k1_select = k1_old[:int(round(pSelfCon*len(X_diag), 0))]
        pairs_select.extend(np.asarray(pairs_diag_old)[k1_select])       
    else:
        N_new = int(round(pSelfCon*len(X_diag), 0)) - N_diag_old
        
        pairs_select.extend(pairs_diag_old)
        # print(pairs_diag)
        # print(type(pairs_diag))
        # print(pairs_diag_old)
        # print(type(pairs_diag_old))
        missing_pairs_id = set(pairs_diag).difference(set(pairs_diag_old))
        k1_new = np.arange(N_new)
        np.random.shuffle(k1_new)
        # print('N new', N_new)
        # print(missing_pairs_id)
        # input('KKK')
        pairs_select.extend(np.asarray(np.asarray(list(missing_pairs_id))[k1_new]))
        
    # print(pairs_select)
    
    idc = Nsyn + np.arange(len(pairs_select)) + 1
    
    pairs_select = np.vstack(pairs_select)
    # print(pairs_select.shape)
    X_out = np.zeros(X.shape)
    X_out[pairs_select[:, 0], pairs_select[:, 1]] = idc

    return X_out, idc
    
    


def find_neigh_recur(X):
    
    Ntarget = X.shape[0]
    Nsource = X.shape[1]
    
    X_local = X.copy()
    
    for i in range(Ntarget):
        X_local[i,i] = 0
    
    X_recur = X_local + X_local.transpose()
    target_in_X, source_in_X = np.where(X_recur>0)
       
    pairs = []
    self_list = [(i, i) for i in range(Nsource)]
    
    for i in range(Ntarget):
        pairs.extend(list(product(source_in_X[target_in_X==i], source_in_X[target_in_X==i])))
        
    pairs_unique, pairs_neighs = np.unique(np.vstack(pairs), axis=0, return_counts=True)
    neigh_mat = np.zeros((Ntarget, Nsource))
    
    for i in range(pairs_unique.shape[0]):
        neigh_mat[pairs_unique[i, 0], pairs_unique[i, 1]]= pairs_neighs[i]
        
    for i in range(Nsource):
        neigh_mat[i,i] = -1
    
    return  neigh_mat
    


    
def p_calc_recur(X, N_neigh, pCon, pSelfCon, slope):
    
    Ntarget = X.shape[0]    
    X_local = X.copy()
    N_neigh_local = np.tril(N_neigh)

    for i in range(Ntarget):
        X_local[i,i] = 0
   
    X_local = np.tril(X_local)
    
    lt = []
    
    for i in range(Ntarget):
        for j in range(i):
            lt.append((i,j))
            
    lt = set(lt)
    pairs_select = set(map(tuple, zip(*np.where(X_local>0))))

    N_neigh_intersect = list(lt.intersection(pairs_select)) 
    N_neigh_ = np.vstack(N_neigh_intersect)
    _, N_neigh_unique = np.unique(N_neigh_, axis=1, return_counts=True)

    N_neigh_selected = N_neigh[N_neigh_[:, 0], N_neigh_[:, 1]]  
    p1 = np.unique(N_neigh[N_neigh>=0])
    p = np.zeros((len(p1), 6))
    p[:, 0] = p1
    
    for i in range(len(p1)):
        pair_act = np.where(N_neigh == p1[i])
        p[i, 1] = len(pair_act[0])
        p[i,2] = np.sum((N_neigh_selected==p1[i]).astype(int))
        p[i,3] = p[i, 2]/p[i, 1]
        
    off = fmin(lambda N0: p_min(N0, p, pCon, slope), np.max(p[:, 0]))

    N1 = np.floor(min(np.max(p[:, 0]), off+1/slope/pCon ))
    N0 = max(np.min(p[:, 0]), np.ceil(off))
    ind = np.isin(p[:, 0], np.arange(N0, N1+1))
    p[ind, 4] = pCon*slope*(p[ind,0]-off)
    ind_type = ind.astype(int)
    
    if len(ind_type) and len(np.where(ind_type)[0]):
        p[np.where(ind_type)[0][-1]+1:, 4] = 1
    else:
        p[:, 4] = 1
    
    N_original = X.shape[0] **2 * pCon - X.shape[0] * pSelfCon
    
    N_act = p[:, 1] @ p[:, 4]
    p[:, 5] = p[:, 4]*N_original/N_act
    
    return p
    
    
def p_min(N0, p, pCon, slope):

    N1 = np.round(min(max(p[:, 0]), N0 + 1/slope/pCon, 0))
    
    N_2 = p[np.where((p[:, 0]>N0)&(p[:, 0]<=N1)), 0]
    pN_2 = p[np.where((p[:, 0]>N0)&(p[:, 0]<=N1)), 1]/np.sum(p[:, 1])
    pN_3 = p[np.where(p[:, 0]>N1), 1]/np.sum(p[:, 1])
    
    return abs(np.sum(pN_2 * slope * (N_2-N0)) + np.sum(pN_3) - 1)

def common_neigh_analysis(X_conn, X_neigh):
    
    X_neigh_ = X_neigh + 1
    X_id = (X_conn>0).astype(int) - (X_conn==0).astype(int)
    X_id = X_neigh_ * X_id - 1
    
    unique1a, unique1b = np.unique(X_id, return_counts=True)
    unique2a, unique2b = np.unique(X_neigh, return_counts=True)

    valid1 = np.where(unique1a>=0)[0].astype(int)
    valid2 = np.where(unique2a>=0)[0].astype(int)
    
    Nneigh1 = unique1a[valid1].astype(int)
    count1 = unique1b[valid1]
    
    Nneigh2 = unique2a[valid2].astype(int)
    count2 = unique2b[valid2]

    Nmax = max(np.max(Nneigh1), np.max(Nneigh2))
    
    Neigharray1 = np.zeros(Nmax+1)
    Neigharray2 = np.zeros(Nmax+1)
    
    Neigharray1[Nneigh1] = count1
    Neigharray2[Nneigh2] = count2
    
    Neigharray = Neigharray1/Neigharray2
    
    return np.arange(Nmax+1), Neigharray
    
def reciprocate_analysis(X):
    
    Ntotal = np.sum((X>0).astype(int))
    Nrec = np.sum((X>0).astype(int) * (X.transpose()>0).astype(int))
    
    return Nrec/Ntotal

def pCon_analysis(X):
    
    return np.sum((X>0).astype(int))/(X.shape[0]*X.shape[1])

if __name__ == '__main__':
    Nsyn=0
    Npop = 200
    pCon = 0.1
    pSelfCon = 0.5
    
    XA, idc = SetCon_CommonNeighbour_Recur(Nsyn, Npop, pCon, pSelfCon)
    
    Npop = 200
    
    XB, idc = SetCon_CommonNeighbour_Recur(Nsyn, Npop, pCon, pSelfCon)
    
 
    
    
    
    
    
