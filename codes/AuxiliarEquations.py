from scipy.integrate import quadrature
import numpy as np
from scipy.signal import argrelextrema as argex


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
    

