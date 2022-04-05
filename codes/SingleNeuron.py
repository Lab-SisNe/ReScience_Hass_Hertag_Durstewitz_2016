from brian2 import *
from AuxiliarEquations import *
from numpy.matlib import repmat
from itertools import product
import os
from scipy.signal import periodogram
from scipy.integrate import quadrature
from scipy.stats import ttest_ind as ttest
from scipy.stats import mannwhitneyu as mwtest


class SingleNeuron():

    def __init__(self, params, I=0, V0=False, w0=0, method='rk4', time_step=0.05):
        
        self.net = Network()
        self.time_step = time_step

        self.params = params
        self.I = I
           
        membane_eq = """
        I: amp
 
        I_ref: amp
        I_exp = g_L * delta_T * exp((V - V_T)/delta_T): amp
        w_V = I + I_exp -g_L * (V - E_L): amp
        
        D0 = (C/g_L) * w_V:  coulomb
        dD0 = C *(exp((V - V_T)/delta_T)-1): farad
        
        dV = int(I >= I_ref) * int(t - last_spike < 5 * ms) * (-g_L/C)*(V - V_dep) + (1 - int(I >= I_ref) * int(t - last_spike < 5 * ms)) * (I + I_exp - g_L * (V - E_L) - w)/C: volt/second
        dV/dt = dV: volt
        dw/dt = int(w > w_V - D0/tau_w) * int(w < w_V + D0/tau_w) * int(V <= V_T) * int(I < I_ref) * -(g_L * (1 - exp((V - V_T)/delta_T)) + dD0/tau_w)*dV: amp
        
        last_spike: second
        C: farad
        g_L: siemens
        E_L: volt
        delta_T:volt
        V_T: volt
        a: siemens
        tau_w: second
        V_up: volt
        V_r: volt
        b: amp
        tau_m: second
        
        V_dep:volt
        """
    
        membane_treshold = "V > V_up"
        membrane_reset = "V = V_r; w += b; last_spike=t"
        membrane_event ='w > w_V - D0/tau_w and w < w_V + D0/tau_w and V <= V_T'

        self.group = NeuronGroup(1, model=membane_eq, threshold=membane_treshold, reset=membrane_reset,
                                 events={'w_crossing': membrane_event}, method=method, refractory=5*ms, dt=self.time_step*ms)
        membrane_eventmonitor = EventMonitor(self.group, "w_crossing", variables=["V", "w"])
        self.group.run_on_event("w_crossing", "w=w_V - D0/tau_w")
        
        
        if type(V0) is int or type(V0) is float:
            _1 = V0
        else:
            _1 = params[2]
        
        if type(w0) is int or type(w0) is float:
            _2 = w0
        else:
            _2 = 0
        
        
        self.group.V = _1*mV
        self.group.w = _2*pA
        
        self.group.I = I*pA

        
        self.group.C = params[0]*pF
        self.group.g_L = params[1]*nS
        self.group.E_L = params[2]*mV
        self.group.delta_T = params[3]*mV
        self.group.V_up = params[4]*mV
        self.group.tau_w = params[5]*ms
        self.group.a = params[6]*nS
        self.group.b = params[7]*pA
        self.group.V_r = params[8]*mV
        self.group.V_T = params[9]*mV
        self.group.I_ref = params[10]*pA
        self.group.V_dep = params[11]*mV
        
        self.group.tau_m = self.group.C/self.group.g_L


        self.I_SN = (self.group.g_L * (self.group.V_T - self.group.E_L - self.group.delta_T))[0]/pA
    
      
        self.spikemonitor = SpikeMonitor(self.group)
        self.neuronmonitor = StateMonitor(self.group, True, True)

        self.net.add(self.group, self.spikemonitor, self.neuronmonitor,
                     membrane_eventmonitor)
      
    def run(self, t):

   
        self.net.run(t, report='text', report_period=10*second)
        
    def plot_V_t(self):
        
        plot(self.neuronmonitor.t/ms, self.neuronmonitor.V[0]/mV)
        tspikes = self.spikemonitor.t/ms
        for t in tspikes:
            vlines(t, self.params[4], 20)
        show()
        
    def plot_w_t(self):
        
        plot(self.neuronmonitor.t/ms, self.neuronmonitor.w[0]/pA)
        tspikes = self.spikemonitor.t/ms
        show()
        
    def phase_portrait(self,I=False, V_range=False, w_range=False):
        
        if type(w_range) is bool:
            w_range = [-1, np.max(self.neuronmonitor.w[0]/pA)+10]
            
        if type(V_range) is bool:
            V_range = [self.params[8]-10, self.params[4] + 10, 0.1]
            
        if type(I) is float or type(I) is int:
            I_ = I
        else:
            I_ = self.I
        V_list = np.arange(V_range[0], V_range[1]+V_range[2], V_range[2])
        
        V_null = np.asarray([w_V(I_, V, self.params) for V in V_list])
        
        plot(V_list, V_null, color='black', label='V null')
        hlines(0, V_range[0], V_range[1], color='black', label='w null')
        vlines(self.params[8], w_range[0], w_range[1], linestyle='--', color='purple', label='V_r')
        vlines(self.params[9], w_range[0], w_range[1], linestyle='--', color='orange', label='V_T')
        vlines(self.params[4], w_range[0], w_range[1], linestyle='--', color='red', label='V_up')
        xlim(V_range[0], V_range[1])
        ylim(w_range[0], w_range[1])
        
        _0 = 0
        for t in self.spikemonitor.t/ms:
            _1 = np.where(self.neuronmonitor.t/ms==t)[0][0]
            plot(self.neuronmonitor.V[0][_0:_1]/mV, self.neuronmonitor.w[0][_0:_1]/pA,  color='blue')
            _0 = _1+1
        else:
            plot(self.neuronmonitor.V[0][_0:]/mV, self.neuronmonitor.w[0][_0:]/pA, label='simulation', color='blue')
            
        
        # plot(self.neuronmonitor.V[0]/mV, self.neuronmonitor.w[0]/pA, label='simulation')
        legend()
        show()
        
        
        

    
if __name__ == '__main__':
    
    # C: 0
    # g_L: 1
    # E_L: 2
    # delta_T: 3
    # V_up: 4
    # tau_w: 5
    # a: 6
    # b:7
    #V_r: 8
    # V_T: 9
    # I_ref = 10
    # V_dep = 11
    
        
    I1 = 56
    I2 = 37
    
    
    print('I1:', I1)
    print('I2:', I2)
    print('-'*40)
    
    print('Not spiking')
    
    params1 = [0 for i in range(12)]
    
    params1[0] = 242.635
    params1[1] = 7.369
    params1[2] = -80.251
    params1[3] = 24.709
    params1[4] = -42.869
    params1[5] = 102.958
    params1[6] = 0
    params1[7] = 6.039
    params1[8] = -67.339
    params1[9] = -48.045
    params1[10] = 1212
    params1[11] = -67.339

    I=I1
    
    neuron1 = SingleNeuron(params1, I=I)
    
    # neuron1.run(5000*ms)
    # neuron1.phase_portrait()
    # neuron1.plot_V_t()
    # neuron1.plot_w_t()
    print('Rheobase:', neuron1.I_SN)
    print('Transient f:', transient_f_I(I, params1))
    print('Steady f:', steady_f_I(I, params1))
    
    print('-'*40)
    
    print('New Not spiking')
    
    params4 = [0 for i in range(12)]
    
    params4[0] = 242.635
    params4[1] = 7.369
    params4[2] = -80.251
    params4[3] = 24.709
    params4[4] = -42.869
    params4[5] = 102.958
    params4[6] = 0
    params4[7] = 6.039
    params4[8] = -72.722
    params4[9] = -48.045
    params4[10] = 1212
    params4[11] = -67.339

    
    neuron4 = SingleNeuron(params4, I=I)
    
    # neuron4.run(5000*ms)
    # neuron4.phase_portrait()
    # neuron4.plot_V_t()
    # neuron4.plot_w_t()
    print('Rheobase:', neuron4.I_SN)
    print('Transient f:', transient_f_I(I, params4))
    print('Steady f:', steady_f_I(I, params4))
    
    print('-'*40)
    print('Spiking')
    
    params2 = [0 for i in range(12)]
    
    params2[0] = 247.943
    params2[1] = 7.456
    params2[2] = -80.634
    params2[3] = 23.163
    params2[4] = -48.416
    params2[5] = 83.243
    params2[6] = 0
    params2[7] = 7.526
    params2[8] = -72.722
    params2[9] = -52.574
    params2[10] = 1212
    params2[11] = -72.722
    
    
    I=I2
    
    neuron2 = SingleNeuron(params2, I=I)
    
    # neuron2.run(5000*ms)
    # neuron2.phase_portrait()
    # neuron2.plot_V_t()
    # neuron2.plot_w_t()
    print('Rheobase:', neuron2.I_SN)
    print('Transient f:', transient_f_I(I, params2))
    print('Steady f:', steady_f_I(I, params2))
    
    
    print('-'*40)
    print('new spiking')
    params3 = [0 for i in range(12)]
    
    params3[0] = 247.943
    params3[1] = 7.456
    params3[2] = -80.634
    params3[3] = 23.163
    params3[4] = -48.416
    params3[5] = 83.243
    params3[6] = 0
    params3[7] = 7.526
    params3[8] = -67.339
    params3[9] = -52.574
    params3[10] = 1212
    params3[11] = -72.722
    
    
    
    
    neuron3 = SingleNeuron(params3, I=I)
    
    # neuron3.run(5000*ms)
    # neuron3.phase_portrait()
    # neuron3.plot_V_t()
    # neuron3.plot_w_t()
    print('Rheobase:', neuron3.I_SN)
    print('Transient f:', transient_f_I(I, params3))
    print('Steady f:', steady_f_I(I, params3))
    
    