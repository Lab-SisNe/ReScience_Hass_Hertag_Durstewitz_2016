from scipy.integrate import quadrature
import numpy as np
from scipy.signal import argrelextrema as argex
from matplotlib import pyplot as plt
from scipy.signal import hilbert, butter, lfilter, periodogram, welch
from scipy.ndimage import gaussian_filter1d as gf1d
import antropy as ant


def w_V(I, V, params):
    
    C, g_L, E_L, delta_T, V_up, tau_w, a, b, V_r, V_T, I_ref, V_dep = params
    tau_m = C/g_L

    return - g_L * (V - E_L) + g_L * delta_T * np.exp((V - V_T)/delta_T) + I

def e_i(I, V, params):
    
    C, g_L, E_L, delta_T, V_up, tau_w, a, b, V_r, V_T, I_ref, V_dep = params
    tau_m = C/g_L

    return w_V(I, V, params) * (1- tau_m/tau_w)

def e_s(I, V, params):
    
    C, g_L, E_L, delta_T, V_up, tau_w, a, b, V_r, V_T, I_ref, V_dep = params
    tau_m = C/g_L

    return w_V(I, V, params) * (1+ tau_m/tau_w)

def w_r(I, params):
    
    C, g_L, E_L, delta_T, V_up, tau_w, a, b, V_r, V_T, I_ref, V_dep = params
    w_th = e_i(I, V_T, params)
    
    return w_th + b

def V_s(I, params, precision=0.001):
    
    C, g_L, E_L, delta_T, V_up, tau_w, a, b, V_r, V_T, I_ref, V_dep = params
    tau_m = C/g_L
    
    w_r_act = w_r(I, params)
        

    if w_r_act < w_V(I, V_r, params) and w_r_act < e_i(I, V_r, params):
        Vi = V_r
        Vf = V_T

        wi = w_r_act - e_i(I, Vi, params)
        wf = w_r_act - e_i(I, Vf, params)

        Vm = np.NaN

        if wi * wf < 0:
            Vm = (Vi + Vf)/2
            wm = w_r_act - e_i(I, Vm, params)

            while Vf - Vi > precision:
                if wm == 0:
                    break
                elif wm * wi > 0:
                    Vi = Vm
                    wi = wm
                else:
                    Vf = Vm

                Vm = (Vi + Vf)/2
                wm = w_r_act - e_i(I, Vm, params)

    elif w_r_act > e_s(I, V_r, params):
        Vi = -100
        Vf = V_r

        wi = w_r_act - e_s(I, Vi, params)
        wf = w_r_act - e_s(I, Vf, params)

        Vm = np.NaN

        if wi * wf < 0:
            Vm = (Vi + Vf)/2
            wm = w_r_act - e_s(I, Vm, params)

            while Vf - Vi > precision:
                if wm == 0:
                    break
                elif wm * wi > 0:
                    Vi = Vm
                    wi = wm
                else:
                    Vf = Vm

                Vm = (Vi + Vf)/2
                wm = w_r_act - e_s(I, Vm, params)
    else:
        Vm = V_r

    return Vm

def V_s_tr(I, params, precision=0.001):
    
    C, g_L, E_L, delta_T, V_up, tau_w, a, b, V_r, V_T, I_ref, V_dep = params
    tau_m = C/g_L
    
    w_r = b

    if w_r < w_V(I, V_r, params) and w_r < e_i(I, V_r, params):
        Vi = V_r
        Vf = V_T

        wi = w_r - e_i(I, Vi, params)
        wf = w_r - e_i(I, Vf, params)

        Vm = np.NaN

        if wi * wf < 0:
            Vm = (Vi + Vf)/2
            wm = w_r - e_i(I, Vm, params)

            while Vf - Vi > precision:

                if wm == 0:
                    break
                elif wm * wi > 0:
                    Vi = Vm
                    wi = wm
                else:
                    Vf = Vm

                Vm = (Vi + Vf)/2
                wm = w_r - e_i(I, Vm, params)

    elif w_r > e_s(I, V_r, params):
        Vi = -100
        Vf = V_r

        wi = w_r - e_s(I, Vi, params)
        wf = w_r - e_s(I, Vf, params)

        Vm = np.NaN

        if wi * wf < 0:
            Vm = (Vi + Vf)/2
            wm = w_r - e_s(I, Vm, params)

            while Vf - Vi > precision:
                if wm == 0:
                    break
                elif wm * wi > 0:
                    Vi = Vm
                    wi = wm
                else:
                    Vf = Vm

                Vm = (Vi + Vf)/2
                wm = w_r - e_s(I, Vm, params)
    else:
        Vm = V_r

    return Vm

def steady_f_I(I, params):
             
    C, g_L, E_L, delta_T, V_up, tau_w, a, b, V_r, V_T, I_ref, V_dep = params
    tau_m = C/g_L
    
    I_SN = g_L * (V_T - E_L - delta_T)
    
    if I <= I_SN:
        return 0
    
    w_r_act = w_r(I, params)
    V_s_act = V_s(I, params)

    def integrand1(V):
        return C/(w_V(I, V, params) - w_r_act)
        
    def integrand2(V):
        return C * tau_w/(tau_m * w_V(I, V, params))

    def integrand3(V):
        return C/(w_V(I, V, params) - w_r_act + b)

    if b > 0:
        t1 = quadrature(integrand1, V_r, V_s_act)[0]
        t2 = quadrature(integrand2, V_s_act, V_T)[0]

        t3 = quadrature(integrand3, V_T, V_up)[0]
        f_ss = 1/(t1 + t2 + t3)

        return f_ss*1000

    elif b == 0:
        t1 = quadrature(integrand1, V_r, V_up)[0]
        f_ss = 1/t1

        return f_ss*1000

def transient_f_I(I, params):

    C, g_L, E_L, delta_T, V_up, tau_w, a, b, V_r, V_T, I_ref, V_dep = params
    tau_m = C/g_L
    
    I_SN = g_L * (V_T - E_L - delta_T)
    
    if I <= I_SN:
        return 0

    def integrand(V):
        return C/(w_V(I, V, params) - b)

    if b <= w_V(I, V_T, params) - tau_m * w_V(I, V_T, params)/tau_w:
        t0 = quadrature(integrand, V_r, V_up)[0]
        f0 = 1/t0

        return f0*1000

    else:
        w_r_act = w_r(I, params)
        V_s = V_s_tr(I, params)
       
        def integrand1(V):
            return C/(w_V(I, V, params) - b)

        def integrand2(V):
            return C * tau_w/(tau_m * w_V(I, V, params))

        def integrand3(V):
            return C/(w_V(I, V, params) - w_r_act + b)

        if b > 0:
            t1 = quadrature(integrand1, V_r, V_s)[0]
            t2 = quadrature(integrand2, V_s, V_T)[0]
          
            t3 = quadrature(integrand3, V_T, V_up)[0]
            
            f_ss = 1/(t1 + t2 + t3)
           

            return f_ss * 1000

        elif b == 0:
            t1 = quadrature(integrand1, V_r, V_up)[0]
            f_ss = 1/t1

            return f_ss *1000
    

def get_min(arr, V_r):
    
    minpos = argex(arr, np.less_equal)
    mins = arr[minpos]
    if np.min(mins) > V_r:
        newmins = mins[mins>V_r]

        arrmin = np.min(newmins)
    else:
        arrmin = np.max(mins)

        
    return arrmin
    
def extract_spikes(arr, V_up, V_r, V_T, t0=0, timestep=0.05):
    

    idc0 = int(t0//timestep + 1)        
    arr = arr[idc0:]   
    
    arrdiff = np.diff(np.abs(arr))
    
    arrmin = max(V_r, np.min(arr))
    ampl = 1 


    sppos = np.where(arrdiff >= ampl)[0][::-1]

    arrmin = get_min(arr, V_r)   

    for i in sppos:
     
        j=i+1
        last = -1000
        spikelist = []
        while arr[j] < arrmin:
            spikelist.append(j)
            if arr[j] - last < 0:
                break
            elif j == len(arr) - 1:
                break
            else:
                last = arr[j]
                j += 1
        
        j = i
        while arr[j] > V_T:
            spikelist.append(j)
            if arr[j] - last > 0:
                break
            elif j==0:
                break
            else:
                last = arr[j]
                j -= 1
        arr = np.delete(arr, spikelist)
        

    return arr       
    
def phase1(sp_arr, time_arr, t0, t1):
    
    # first_t_arr = np.arange(t0, t1, dt)
    phase_arr = np.ones(len(time_arr)) * (-1)

    for p in range(len(sp_arr)-1):
        
        incycle = np.where((time_arr>=sp_arr[p]) & (time_arr <sp_arr[p+1]))[0]

        phase_arr[incycle] = np.linspace(p*2 *np.pi,(p+1)* 2*np.pi, len(incycle)) 

    t_bool = np.where((time_arr>=t0)&(time_arr<t1))[0]
    time_arr = time_arr[t_bool]
    phase_arr = phase_arr[t_bool]
    
    return time_arr, phase_arr


def PLV_from_phase_st(ph1, ph2):
    
    bool1 = np.where(ph1>-1)[0]
    ph1 = ph1[bool1]
    ph2 = ph2[bool1]
    bool2 = np.where(ph2>-1)[0]
    ph1 = ph1[bool2]
    ph2 = ph2[bool2]

    diffph = ph2 - ph1 
    exp_sum = np.sum(np.exp(1j * diffph))
    
    return np.absolute(exp_sum)/len(diffph)

def PLV_st(sp_list, time_arr, t0, t1):
    
    phase_list = []
    for i in range(len(sp_list)):
        
        phase_list.append(phase1(sp_list[i], time_arr, t0, t1)[1])
        
    # input()
    pop = []
    phased = []
    for i in range(len(phase_list)):
        for k in range(0, i):
            _ = PLV_from_phase_st(phase_list[i], phase_list[k])
            if not np.isnan(_):
                pop.append(_)
                
    return np.mean(pop)

def phase_signal(sign_arr, time_arr, t0, t1):
    
    bool_arr = np.where((time_arr>=t0) & (time_arr<t1))[0]
    t_arr = time_arr[bool_arr] - t0
    signal = sign_arr[bool_arr]
    
    analytic_signal = hilbert(signal)
    instantaneous_phase = np.unwrap(np.angle(analytic_signal))
    
    return t_arr, instantaneous_phase             

def PLV_from_phase_diff(theta1, theta2):
    complex_phase_diff = np.exp(complex(0,1)*(theta1 - theta2))
    plv = np.abs(np.sum(complex_phase_diff))/len(theta1)
    return plv

def PLV_signal(sgn_list, time_list, t0, t1):
    
    phase_list = []
    for i in range(len(sgn_list)):
        phase_list.append(phase_signal(sgn_list[i], time_list, t0, t1)[1])
    
    pop = []
    for i in range(len(phase_list)):
        for k in range(0, i):
            
            pop.append(PLV_from_phase_diff(phase_list[i],phase_list[k]))
  
    return np.mean(pop)

def PLV_signal_filtered(sgn_list, time_list, t0, t1, lowcut, highcut, fs=20000, order=1):
    filtered_sgn_list = []
    
    for i in range(len(sgn_list)):
        filtered_sgn_list.append(butter_bandpass_filter(sgn_list[i], lowcut, highcut, fs, order))
    
    return PLV_signal(filtered_sgn_list, time_list, t0, t1)

def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filter(data, lowcut, highcut, fs=20000, order=1):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def arrange_spikes(sp_list, idc_list, idc):
    
    idc_bool = np.isin(idc_list, idc)
    idc_list = idc_list[idc_bool]
    sp_list = sp_list[idc_bool]
    
    new_list = []
    
    for i in idc:
        newsp = sp_list[idc_list==i]
        if len(newsp)>0:
            new_list.append(newsp)
    
    
    return new_list
    
def pop_rate(sp_list, t0, t1, dt, sigma=False):
    
    Ntbins = int(round((t1-t0)/dt, 0))
    Nsparrays = len(sp_list)
    
    sp_array = np.zeros((Nsparrays, Ntbins))
    t = np.arange(t0, t1, dt)
    
    for i in range(Nsparrays):
        sparr =np.asarray(sp_list[i])
        spbool = np.where((sparr>=t0) & (sparr<t1))[0]
        sparr = sparr[spbool] - t0
        idc_arr = (sparr//dt).astype(int)

        sp_array[i, idc_arr] = 1
    
    rate_arr = np.mean(sp_array, axis=0)*1000/dt
    
    if type(sigma) is not bool:
        rate_arr = gf1d(rate_arr, sigma=sigma)
    
    return t, rate_arr




def xlogx(x, base=2):
    """Returns x log_b x if x is positive, 0 if x == 0, and np.nan
    otherwise. This handles the case when the power spectrum density
    takes any zero value.
    """
    x = np.asarray(x)
    xlogx = np.zeros(x.shape)
    xlogx[x < 0] = np.nan
    valid = x > 0
    xlogx[valid] = x[valid] * np.log(x[valid]) / np.log(base)
    return xlogx


def SE_signal(signal_array, t0, t1, dt, frequency=False, return_periodogram=False):
    sf = 1000/dt
    time_arr = np.arange(t0, t1, dt)
    time_bool = np.where((time_arr>=t0) & (time_arr<t1))[0]
    signal_array = signal_array[time_bool]
    
    fq, psd = periodogram(signal_array, sf)
    
    if type(frequency) is not bool:
        f0, f1 = frequency
        f_bool = np.where((fq>= f0)&(fq<f1))[0]
        fq = fq[f_bool]
        psd = psd[f_bool]
    
    psd_norm = psd / psd.sum()
    se = -xlogx(psd_norm).sum()

    se /= np.log2(len(psd_norm))
    if return_periodogram:
        return se, fq, psd
    else:
        return se
    

def SE_signal_group1(signal_matrix, t0, t1, dt, not_matrix=False, return_periodogram=False):
    """ Returns network spectral entropy as the spectral entropy of the mean SPD"""
    
    
    
    if not_matrix:
        mat = np.zeros((len(signal_matrix), len(signal_matrix[0])))
        for i in range(len(signal_matrix)):
            mat[i] = signal_matrix[i]
        signal_matrix = mat
    
    sf = 1000/dt
    time_arr = np.arange(t0, t1, dt)
    time_bool = np.where((time_arr>=t0) & (time_arr<t1))[0]
    signal_matrix = signal_matrix[:, time_bool]
    
    period_list = np.zeros((len(signal_matrix), int(len(signal_matrix[0])//2) + 1))
    for x in range(len(signal_matrix)):
        fq, period_list[x] = periodogram(signal_matrix[x], sf)
    psd = np.mean(period_list, axis=0)
    
    psd_norm = psd / psd.sum()
    se = -xlogx(psd_norm).sum()

    se /= np.log2(len(psd_norm))
    if return_periodogram:
        return se, fq, psd
    else:
        return se

def SE_signal_group2(signal_matrix, t0, t1, dt, not_matrix=False):
    """ Returns network spectral entropy as the mean of individual SPD spectral entropy"""
    
    if not_matrix:
        mat = np.zeros((len(signal_matrix), len(signal_matrix[0])))
        for i in range(len(signal_matrix)):
            mat[i] = signal_matrix[i]
        signal_matrix = mat
    
    time_arr = np.arange(t0, t1, dt)
    time_bool = np.where((time_arr>=t0) & (time_arr<t1))[0]
    signal_matrix = signal_matrix[:, time_bool]
    
    SE = []
    fs = 1000/dt
    for i in range(signal_matrix.shape[0]):
        SE.append(ant.spectral_entropy(signal_matrix[i], fs, normalize=True))
        
    SE = np.asarray(SE)
    
    return np.mean(SE)

def set_binary_vector(spike_times, t0, t1, dt):
    
    spike_times = np.asarray(spike_times)
 
    t_bool = np.where((spike_times>=t0)&(spike_times<t1))[0] 
    spike_times = spike_times[t_bool]-t0
    t_bins = (spike_times/dt).astype(int)
    
    t_vector= np.zeros(int((t1-t0)/dt))
    t_vector[t_bins] = 1
    
    return t_vector

def p_vector(t_vector):
    
    return np.sum(t_vector)/len(t_vector)

def p_joint_vector(t_vec1, t_vec2, lag):
    
    if lag>0:
        vector1 = t_vec1[:-lag]
        vector2 = t_vec2[lag:]
        
        return (vector1 @ vector2)/len(vector1)
    
    elif lag<0:
        vector1 = t_vec1[-lag:]
        vector2 = t_vec2[:lag]
        
        return (vector1 @ vector2)/len(vector1)
    
    else:
        return (t_vec1 @ t_vec2)/len(t_vec1)
    
def Pearson(t_vec1, t_vec2, lag):
    t_vec1 = np.asarray(t_vec1)
    t_vec2 = np.asarray(t_vec2)
    
    if lag >0:
        vector1 = t_vec1[:-lag]
        vector2 = t_vec2[lag:]
        
    elif lag <0:
        vector1 = t_vec1[-lag:]
        vector2 = t_vec2[:lag]
        
    else:
        vector1 = t_vec1
        vector2 = t_vec2
    
    px = p_vector(vector1)
    py = p_vector(vector2)
    
    pxy = p_joint_vector(t_vec1, t_vec2, lag)
    
    if px==0  or py ==0:
        return 0
    else:
        return (pxy - px * py)/np.sqrt(px*(1-px)*py*(1-py))

def Pearson_correlation(spike_times1, spike_times2, t0, t1, t_bin, lag0, lag1):
    
    
    t_vec1 = set_binary_vector(spike_times1, t0, t1, t_bin)
    t_vec2 = set_binary_vector(spike_times2, t0, t1, t_bin)
    
    lag_array = np.arange(lag0, lag1, t_bin)
    lag_array_dt = (lag_array/t_bin).astype(int)
    
    corr = []
    # print(t_vec1)
    # print(t_vec2)
    for lag in lag_array_dt:
        corr.append(Pearson(t_vec1, t_vec2, lag))
    
    return lag_array, np.asarray(corr)

def Pearson_correlation_group(spike_group, t0, t1, t_bin, lag0, lag1):
    
    Nlag = int((lag1 -lag0)/t_bin)
    Ngroup = len(spike_group)
    
    corr_group = np.zeros((int(Ngroup*(Ngroup-1)/2), Nlag))
    k = 0
    for i in range(Ngroup):
        for j in range(0, i):
            
            spike_times1 = spike_group[i]
            spike_times2 = spike_group[j]
            
            t_vec1 = set_binary_vector(spike_times1, t0, t1, t_bin)
            t_vec2 = set_binary_vector(spike_times2, t0, t1, t_bin)
    
            corr_group[k] = Pearson_correlation(spike_times1, spike_times2, t0, t1, t_bin, lag0, lag1)[1]
            
            k+= 1
    
    lag_array = np.arange(lag0, lag1, t_bin)
    return lag_array, np.mean(corr_group, axis=0)
