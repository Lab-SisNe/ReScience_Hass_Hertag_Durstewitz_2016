input_ = input
from brian2 import *
from AuxiliarEquations import *
from numpy.matlib import repmat
from itertools import product
import os
from scipy.signal import periodogram
from scipy.integrate import quadrature
from scipy.stats import ttest_ind as ttest
from scipy.stats import mannwhitneyu as mwtest
from scipy.stats import chi2_contingency as chi2

class CortexNetwork():

    def __init__(self, NeuPar,V0, STypPar,SynPar,SPMtx, group_distr,
                 constant_current,fluctuating_current,
                 scales2,
                 method, time_step,          
                 simseed):
        
        seed(simseed)
        
        self.net = Network()
        self.time_step = time_step

        self.group_distr = group_distr
        self.NeuPar = NeuPar
        self.SynPar = SynPar
        
        sourceAMPA_gmaxscale, sourceGABA_gmaxscale, sourceNMDA_gmaxscale, targetAMPA_gmaxscale, targetGABA_gmaxscale, targetNMDA_gmaxscale = scales2
        
        tau_on_AMPA = STypPar[1, 0]*ms
        tau_on_GABA = STypPar[1, 1]*ms
        tau_on_NMDA = STypPar[1, 2]*ms
        
        tau_off_AMPA = STypPar[2, 0]*ms
        tau_off_GABA = STypPar[2, 1]*ms
        tau_off_NMDA =  STypPar[2, 2]*ms
        
        E_AMPA = STypPar[3, 0] * mV
        E_GABA = STypPar[3, 1] * mV
        E_NMDA = STypPar[3, 2] * mV

        Mg_fac = STypPar[5, 2]
        Mg_slope = STypPar[6, 2]
        Mg_half = STypPar[7, 2]

        self.Gsyn_AMPA =  STypPar[0, 0] * STypPar[2, 0] * STypPar[1, 0]/(STypPar[2, 0] - STypPar[1, 0])
        self.Gsyn_GABA =  STypPar[0, 1] * STypPar[2, 1] * STypPar[1, 1]/(STypPar[2, 1] - STypPar[1, 1])
        self.Gsyn_NMDA =  STypPar[0, 2] * STypPar[2, 2] * STypPar[1, 2]/(STypPar[2, 2] - STypPar[1, 2])

        if len(fluctuating_current):
          
            start_fluct, stop_fluct, AC = fluctuating_current
            
            nsteps = int(round(stop_fluct/time_step, 0))
            nstart = int(round(start_fluct/time_step, 0))
            
            self.ta = np.zeros((nsteps, NeuPar.shape[1]))
            
            for i in range(len(group_distr[0])):
                for j in range(len(group_distr)):
                    fluct_str = AC[i][j]
                    if fluct_str != '':
                        I_arr = eval("asarray([0 for i in range(nstart)]+[{} for t in linspace({}, {}, {}, endpoint=False)])".format(fluct_str, start_fluct, stop_fluct, nsteps-nstart))
                        cells = np.asarray(group_distr[j][i])
                        subta = repmat(I_arr, len(cells), 1).transpose()
                        self.ta[:, cells] = subta
                        
            eq_fluc = 'ta(t, i)'
            
        else:
            self.ta = []
            eq_fluc = '0'
    
        membane_eq = """
        I_DC: amp
        I_AC = {}*pA: amp
        I_inj = I_DC + I_AC: amp      
        I_syn = I_AMPA + I_NMDA + I_GABA: amp
        I_tot =  I_syn + I_inj: amp
        
        I_EXC = I_AMPA + I_NMDA: amp
        I_EXC_tot = I_EXC + I_inj: amp
        
        I_ref: amp
        I_exp = g_L * delta_T * exp((V - V_T)/delta_T): amp
        w_V = I_tot + I_exp -g_L * (V - E_L): amp
        
        D0 = (C/g_L) * w_V:  coulomb
        dD0 = C *(exp((V - V_T)/delta_T)-1): farad
        
        dV = int(I_tot >= I_ref) * int(t - last_spike < 5 * ms) * (-g_L/C)*(V - V_dep) + (1 - int(I_tot >= I_ref) * int(t - last_spike < 5 * ms)) * (I_tot + I_exp - g_L * (V - E_L) - w)/C: volt/second
        dV/dt = dV: volt
        dw/dt = int(w > w_V - D0/tau_w) * int(w < w_V + D0/tau_w) * int(V <= V_T) * int(I_tot < I_ref) * -(g_L * (1 - exp((V - V_T)/delta_T)) + dD0/tau_w)*dV: amp
        
        
        dg_AMPA_off/dt = - (1/tau_off_AMPA) * g_AMPA_off: siemens
        dg_AMPA_on/dt = - (1/tau_on_AMPA) * g_AMPA_on: siemens
        
                
        dg_GABA_off/dt = -(1/tau_off_GABA)*  g_GABA_off: siemens
        dg_GABA_on/dt = -(1/tau_on_GABA) *  g_GABA_on: siemens 
        
        dg_NMDA_off/dt = -(1/tau_off_NMDA)*  g_NMDA_off: siemens
        dg_NMDA_on/dt = -(1/tau_on_NMDA) *  g_NMDA_on: siemens
        
        g_AMPA = g_AMPA_off - g_AMPA_on: siemens
        g_GABA = g_GABA_off - g_GABA_on: siemens
        g_NMDA = (g_NMDA_off - g_NMDA_on): siemens
        
        I_AMPA = (g_AMPA) * (E_AMPA - V): amp
        Mg = (1/(1 +  {} * exp({} * ( {}*mV - V)/mV))):1
        I_GABA = (g_GABA) *(E_GABA - V): amp
        I_NMDA = (g_NMDA) * (E_NMDA - V)* Mg: amp        
        	
        E_AMPA = {} * volt: volt
        E_GABA = {} * volt: volt
        E_NMDA = {} * volt: volt
        
        tau_on_AMPA = {}*second: second
        tau_on_GABA = {}*second: second
        tau_on_NMDA = {}*second: second
        
        tau_off_AMPA = {}*second: second
        tau_off_GABA = {}*second: second
        tau_off_NMDA =  {}*second: second
    
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
        """.format(eq_fluc, Mg_fac, Mg_slope, Mg_half, E_AMPA, E_GABA, E_NMDA,
        tau_on_AMPA, tau_on_GABA, tau_on_NMDA, tau_off_AMPA, tau_off_GABA, tau_off_NMDA)
    
        membane_treshold = "V > V_up"
        membrane_reset = "V = V_r; w += b"
        membrane_event ='w > w_V - D0/tau_w and w < w_V + D0/tau_w and V <= V_T'

        self.group = NeuronGroup(NeuPar.shape[1], model=membane_eq, threshold=membane_treshold, reset=membrane_reset,
                                 events={'w_crossing': membrane_event}, method=method, refractory=5*ms, dt=self.time_step*ms)
        membrane_eventmonitor = EventMonitor(self.group, "w_crossing", variables=["V", "w"])
        self.group.run_on_event("w_crossing", "w=w_V - D0/tau_w")
        
        
        if len(constant_current):
            for i in range(len(group_distr[0])):
                for j in range(len(group_distr)):
                    _ = np.asarray(group_distr[j][i]).astype(int)                
                    self.group.I_DC[_] = constant_current[i][j]*pA
        
        self.group.last_spike = -100000000*ms
       
        self.group.C = NeuPar[0, :]*pF
        self.group.g_L = NeuPar[1, :]*nS
        self.group.E_L = NeuPar[2, :]*mV
        self.group.delta_T = NeuPar[3, :]*mV
        self.group.V_up = NeuPar[4, :]*mV
        self.group.tau_w = NeuPar[5, :]*ms
        self.group.a = NeuPar[6, :]*nS
        self.group.b = NeuPar[7, :]*pA
        self.group.V_r = NeuPar[8, :]*mV
        self.group.V_T = NeuPar[9, :]*mV
        self.group.I_ref = NeuPar[10, :]*pA
        self.group.V_dep = NeuPar[11, :]*mV
        
        self.group.V = V0[0] * mV
        self.group.w = V0[1] * pA

        self.I_SN = self.group.g_L * (self.group.V_T - self.group.E_L - self.group.delta_T)
    
        print('REPORT: NeuronGroup set\n')
        
        synapse_model = """
        AMPA: 1
        GABA: 1
        NMDA: 1
        u1: 1
        u: 1
        U: 1
        R: 1
        R1: 1
        failure: 1
        failure_rate: 1
        a_syn: 1
        tau_fac: second
        tau_rec: second
        g_max: siemens
        Gsyn: 1
        """
        syn_pathway0 = "failure = int(rand()<failure_rate)"
        syn_pathway1 = "u1 = U + u * (1 - U) * exp(-(t - last_spike_pre)/tau_fac)"
        syn_pathway2 = "R1 = 1 + (R - u * R - 1) * exp(- (t - last_spike_pre)/tau_rec)"
        syn_pathway3 = 'u = u1'
        syn_pathway4 = 'R = R1'
        syn_pathway5 = 'last_spike_pre = t'
        syn_pathway6 = 'a_syn = u * R'
        syn_pathway7a = "g_AMPA_on_post += AMPA * g_max * a_syn * (1 - failure) * Gsyn"
        syn_pathway7b = "g_AMPA_off_post += AMPA * g_max * a_syn * (1 - failure) * Gsyn"
        syn_pathway8a = "g_GABA_on_post += GABA * g_max * a_syn * (1 - failure) * Gsyn"
        syn_pathway8b = "g_GABA_off_post += GABA * g_max * a_syn * (1 - failure) * Gsyn"
        syn_pathway9a = "g_NMDA_on_post += NMDA * g_max * a_syn * (1 - failure) * Gsyn"
        syn_pathway9b = "g_NMDA_off_post += NMDA* g_max * a_syn * (1 - failure) * Gsyn"
        
        syn_dict_pathway = {'p0': syn_pathway0,   
                             'p1': syn_pathway1,
                             'p2': syn_pathway2,
                             'p3': syn_pathway3,
                             'p4': syn_pathway4,
                             'p5': syn_pathway5,
                             'p6': syn_pathway6,
                             'p7a': syn_pathway7a,
                             'p7b': syn_pathway7b,
                             'p8a': syn_pathway8a,
                             'p8b': syn_pathway8b,
                             'p9a': syn_pathway9a,
                             'p9b': syn_pathway9b,
                             }
  
        self.target_arr, self.source_arr = sourcetarget(SPMtx)

        self.synapses = Synapses(self.group, self.group, model=synapse_model, on_pre=syn_dict_pathway, method=method, dt=self.time_step*ms)
        
        if len(self.target_arr):
            self.synapses.connect(i=self.source_arr, j=self.target_arr)
            
            stype = SynPar[0, :]
            use = SynPar[1, :]
            tau_rec = SynPar[2, :]*ms
            tau_fac = SynPar[3, :]*ms
            g_max = SynPar[4, :]*nS
            dtax = SynPar[5, :]*ms
            p_fail = SynPar[6, :]
            R  =  1
            a = use
            
            self.synapses.AMPA = ((3- stype) * (2 - stype)//2)
            self.synapses.GABA = (stype-1) * (3 - stype)
            self.synapses.NMDA = (stype-1)*(stype-2)//2
            
            Gsyn_AMPA_arr = self.synapses.AMPA * self.Gsyn_AMPA
            Gsyn_GABA_arr = self.synapses.GABA * self.Gsyn_GABA
            Gsyn_NMDA_arr = self.synapses.NMDA * self.Gsyn_NMDA
            
            self.synapses.Gsyn = Gsyn_AMPA_arr + Gsyn_GABA_arr + Gsyn_NMDA_arr
           
            self.synapses.U = use
            self.synapses.u = use
            self.synapses.u1 = use
            self.synapses.R = R
            self.synapses.R1 = R
            self.synapses.a_syn = a
            
            self.synapses.tau_rec = tau_rec
            self.synapses.tau_fac = tau_fac
            self.synapses.g_max = g_max
            self.synapses.failure_rate = p_fail
            
            
            sourceAMPA_gmaxscale, sourceGABA_gmaxscale, sourceNMDA_gmaxscale, targetAMPA_gmaxscale, targetGABA_gmaxscale, targetNMDA_gmaxscale = scales2

            for group, stripe, scale in sourceAMPA_gmaxscale:
                AMPAsyn = self.group_name_to_AMPAsynapses_idc([['all', 'all', group, stripe]])
                self.synapses.g_max[AMPAsyn] *= scale
            
            for group, stripe, scale in sourceGABA_gmaxscale:
                GABAsyn = self.group_name_to_GABAsynapses_idc([['all', 'all', group, stripe]])
                self.synapses.g_max[GABAsyn] *= scale
                
            for group, stripe, scale in sourceNMDA_gmaxscale:
                NMDAsyn = self.group_name_to_NMDAsynapses_idc([['all', 'all', group, stripe]])
                self.synapses.g_max[NMDAsyn] *= scale
                
            for group, stripe, scale in targetAMPA_gmaxscale:
                AMPAsyn = self.group_name_to_AMPAsynapses_idc([[group, stripe, 'all', 'all']])
                a = self.synapses.g_max[AMPAsyn]
                self.synapses.g_max[AMPAsyn] *= scale
    
            
            for group, stripe, scale in targetGABA_gmaxscale:
                GABAsyn = self.group_name_to_GABAsynapses_idc([[group, stripe, 'all', 'all']])
                self.synapses.g_max[GABAsyn] *= scale
                
            for group, stripe, scale in targetNMDA_gmaxscale:
                NMDAsyn = self.group_name_to_NMDAsynapses_idc([[group, stripe, 'all', 'all']])
                a = self.synapses.g_max[NMDAsyn]
                self.synapses.g_max[NMDAsyn] *= scale
                
            
            self.synapses.p7a.delay = dtax
            self.synapses.p7b.delay = dtax
            self.synapses.p8a.delay = dtax
            self.synapses.p8b.delay = dtax
            self.synapses.p9a.delay = dtax
            self.synapses.p9b.delay = dtax

            self.synapses.p0.order = 0
            self.synapses.p1.order = 1
            self.synapses.p2.order = 2
            self.synapses.p3.order = 3
            self.synapses.p4.order = 4
            self.synapses.p5.order = 5
            self.synapses.p6.order = 6
            self.synapses.p7a.order = 7
            self.synapses.p7b.order = 7
            self.synapses.p8a.order = 7
            self.synapses.p8b.order = 7
            self.synapses.p9a.order = 7
            self.synapses.p9b.order = 7
        
        else:
            self.synapses.active = False
            
        self.spikemonitor = SpikeMonitor(self.group)
        self.neuronmonitors = {}
        self.synapsesmonitors = {}
        
        self.net.add(self.group, self.synapses, self.spikemonitor,
                     membrane_eventmonitor)
      
        print('REPORT: Synapses set\n')  
        
        print('REPORT: Network set\n')
        
    def group_name_to_group_idc(self, group):
    
        group_dict = {'PC_L23': [0,], 'IN_L_L23': [1,], 'IN_L_d_L23': [2,], 'IN_CL_L23': [3,], 'IN_CL_AC_L23': [4,], 'IN_CC_L23': [5,], 'IN_F_L23': [6,],
                             'PC_L5': [7,], 'IN_L_L5': [8,], 'IN_L_d_L5': [9,], 'IN_CL_L5': [10,], 'IN_CL_AC_L5':[11,], 'IN_CC_L5': [12,], 'IN_F_L5': [13,],
                             'PC': [0, 7], 'IN': [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13],
                             'IN_L': [1, 8], 'IN_L_d': [2, 9],'IN_CL': [3, 10], 'IN_CL_AC': [4, 11], 'IN_CC': [5, 12], 'IN_F': [6, 13],
                             'IN_L23': [1, 2, 3, 4, 5, 6], 'IN_L5': [8, 9, 10, 11, 12, 13],
                             'L_23': [i for i in range(7)], 'L_5': [i for i in range(7, 14)],
                             'all': [i for i in range(14)],
                             }
        return group_dict[group]
    
    def group_idc_to_group_name(self, idc):
        
        idc_dict = {0: 'PC_L23', 1: 'IN_L_L23', 2: 'IN_L_d_L23', 3: 'IN_CL_L23', 4: 'IN_CL_AC_L23', 5: 'IN_CC_L23', 6: 'IN_F_L23',
                    7: 'PC_L5', 8: 'IN_L_L5', 9: 'IN_L_d_L5', 10: 'IN_CL_L5', 11: 'IN_CL_AC_L5', 12: 'IN_CC_L5', 13: 'IN_F_L5',
                    }
        
        return idc_dict[idc]
    
    def neuron_idc_to_group_name(self, idc):
        
        for group in range(len(self.group_distr)):
            for stripe in range(len(self.group_distr[0])):
                if idc in self.group_distr[group][stripe]:
                    return self.group_idc_to_group_name(group), stripe
       
    def neuron_monitors(self, MonitorsSetup):
                
        for monitor in MonitorsSetup:
            if len(monitor[2]) == 1 and monitor[2][0][0] == 'all' and monitor[2][0][1] == 'all':
                neuron_idc = True
            else:
                neuron_idc = self.group_name_to_neuron_idc(monitor[2])
                
            self.neuronmonitors[monitor[0]] = StateMonitor(self.group, monitor[1], neuron_idc, dt=self.time_step*ms)
       
        self.net.add(*self.neuronmonitors.values())
        
    def synapses_monitors(self, MonitorsSetup):
        
        
        for monitor in MonitorsSetup:         
            syn_idc = self.group_name_to_synapses_idc(monitor[2])
            self.synapsesmonitors[monitor[0]] = StateMonitor(self.group, monitor[1], syn_idc, dt=self.time_step*ms)
       
        self.net.add(*self.synapsesmonitors.values())
            
    def poisson_input(self, PoissonInput):
        
        poisson_source = []
        poisson_target = []
        poisson_type = []
        poisson_gmax = []
        time_poisson = []
        spiking_poisson = []
        
        NPoissonTarget = 0
        NPoissonInp = 0
 
        for inp in PoissonInput:    
            Ninp, Typeinp, Freqinp, gmaxinp, pfailinp, Startinp, Stopinp, targetsinp = inp
            inp_group = np.arange(NPoissonInp, NPoissonInp+Ninp)
            
            for i in range(Ninp):
                Lambinp = 1000/Freqinp
                N_inp = int(round(2*(Stopinp-Startinp)/Lambinp, 0))
                spikepoisson =  Startinp + np.cumsum(np.clip(-Lambinp *np.log(1-np.random.rand(N_inp)), self.time_step, np.NaN))
                spikeact = spikepoisson[spikepoisson < Stopinp]
                        
                time_poisson.extend(spikeact)
                spiking_poisson.extend([NPoissonInp+i for sp in spikeact])
                
                for target in targetsinp:
                    groupinp, stripeinp, pConinp = target
                    
                    target_group = self.group_distr[groupinp][stripeinp]
                    Ntarget = int(round(pConinp*len(target_group), 0))

                    comb_inp = list(product(target_group, inp_group))
                    np.random.shuffle(comb_inp)
                    select_inp = comb_inp[:Ntarget]
                    select_target, select_source = list(map(list, zip(*select_inp)))
                    poisson_source.extend(select_source)
                    poisson_target.extend(select_target)
                    poisson_type.extend([Typeinp for i in range(Ntarget)])
                    poisson_gmax.extend([gmaxinp for i in range(Ntarget)])
                    
            NPoissonInp += Ninp
           
        self.poissoninput = SpikeGeneratorGroup(NPoissonInp, spiking_poisson, time_poisson*ms, dt=self.time_step*ms )
        
        poisson_synapse = """
        AMPA: 1
        GABA: 1
        NMDA: 1
        failure: 1
        failure_rate: 1
        g_max: siemens
        Gsyn_AMPA: 1
        Gsyn_GABA: 1
        Gsyn_NMDA: 1
        """
        poisson_pathway0 = "failure =  int(rand()<failure_rate)"
        poisson_pathway7a = "g_AMPA_on_post += AMPA * g_max * (1 - failure) * Gsyn_AMPA"
        poisson_pathway7b = "g_AMPA_off_post += AMPA * g_max * (1 - failure) * Gsyn_AMPA"
        poisson_pathway8a = "g_GABA_on_post += GABA * g_max *  (1 - failure) * Gsyn_GABA"
        poisson_pathway8b = "g_GABA_off_post += GABA * g_max *  (1 - failure) * Gsyn_GABA"
        poisson_pathway9a = "g_NMDA_on_post += NMDA * g_max *  (1 - failure) * Gsyn_NMDA"
        poisson_pathway9b = "g_NMDA_off_post += NMDA* g_max *  (1 - failure) * Gsyn_NMDA"
        
        poisson_dict = {'p0': poisson_pathway0,
                             'p7a': poisson_pathway7a,
                             'p7b': poisson_pathway7b,
                             'p8a': poisson_pathway8a,
                             'p8b': poisson_pathway8b,
                             'p9a': poisson_pathway9a,
                             'p9b': poisson_pathway9b,
                             }
        
        self.poissonsynapse = Synapses(self.poissoninput, self.group, model=poisson_synapse, on_pre=poisson_dict, dt=self.time_step*ms)
        
        poisson_source = np.asarray(poisson_source).astype(int)
        poisson_target = np.asarray(poisson_target).astype(int)
        
        self.poissonsynapse.connect(i=poisson_source, j=poisson_target)
        
        poisson_type = np.asarray(poisson_type)
        
        self.poissonsynapse.AMPA = 2 - poisson_type
        self.poissonsynapse.GABA = poisson_type - 1
        self.poissonsynapse.NMDA = 2 - poisson_type
        
        self.poissonsynapse.failure_rate = pfailinp
        self.poissonsynapse.g_max = poisson_gmax *nS
        
        self.poissonsynapse.Gsyn_AMPA = self.Gsyn_AMPA
        self.poissonsynapse.Gsyn_GABA = self.Gsyn_GABA
        self.poissonsynapse.Gsyn_NMDA = self.Gsyn_NMDA
        
        self.poissonspikemonitor = SpikeMonitor(self.poissoninput)
    
        self.net.add(self.poissoninput, self.poissonsynapse, self.poissonspikemonitor)

    def regular_input(self, RegularInput):
        
        regular_source = []
        regular_target = []
        regular_type = []

        regular_gmax = []
        time_regular = []
        spiking_regular = []
        
        NregularTarget = 0
        NregularInp = 0
        
        for inp in RegularInput:    
            Ninp, Typeinp, Nspikes, gmaxinp, pfailinp, Startinp, Stopinp, targetsinp = inp
            inp_group = np.arange(NregularInp, NregularInp+Ninp)
            
            if (Stopinp-Startinp)/(Nspikes-1)< self.time_step:
                print('WARNING: time step is larger than spikes intervals')
                
            for i in range(Ninp):
                time_regular.extend(np.linspace(Startinp, Stopinp, Nspikes))
                spiking_regular.extend([NregularInp+i for sp in range(Nspikes)])
            
                for target in targetsinp:
                    groupinp, stripeinp, pConinp = target
                    
                    target_group = self.group_distr[groupinp][stripeinp]
                    Ntarget = int(round(pConinp*len(target_group), 0))
       
                    comb_inp = list(product(target_group, inp_group))
                    np.random.shuffle(comb_inp)
                    select_inp = comb_inp[:Ntarget]
                    select_target, select_source = list(map(list, zip(*select_inp)))
                    regular_source.extend(select_source)
                    regular_target.extend(select_target)
                    regular_type.extend([Typeinp for i in range(Ntarget)])
                    regular_gmax.extend([gmaxinp for i in range(Ntarget)])
                    
            NregularInp += Ninp

        self.regularinput = SpikeGeneratorGroup(NregularInp, spiking_regular, time_regular*ms, dt=self.time_step*ms )
        
        regularmodel = """
        AMPA: 1
        GABA: 1
        NMDA: 1
        failure: 1
        failure_rate: 1
        g_max: siemens
        Gsyn_AMPA: 1
        Gsyn_GABA: 1
        Gsyn_NMDA: 1
        """
        regular_pathway0 = "failure =  int(rand()<failure_rate)"
        regular_pathway7a = "g_AMPA_on_post += AMPA * g_max * (1 - failure) * Gsyn_AMPA"
        regular_pathway7b = "g_AMPA_off_post += AMPA * g_max * (1 - failure) * Gsyn_AMPA"
        regular_pathway8a = "g_GABA_on_post += GABA * g_max *  (1 - failure) * Gsyn_GABA"
        regular_pathway8b = "g_GABA_off_post += GABA * g_max *  (1 - failure) * Gsyn_GABA"
        regular_pathway9a = "g_NMDA_on_post += NMDA * g_max *  (1 - failure) * Gsyn_NMDA"
        regular_pathway9b = "g_NMDA_off_post += NMDA* g_max *  (1 - failure) * Gsyn_NMDA"
        
        regular_dict = {'p0': regular_pathway0,
                             'p7a': regular_pathway7a,
                             'p7b': regular_pathway7b,
                             'p8a': regular_pathway8a,
                             'p8b': regular_pathway8b,
                             'p9a': regular_pathway9a,
                             'p9b': regular_pathway9b,
                             }
        self.regularsynapse = Synapses(self.regularinput, self.group, model=regularmodel, on_pre=regular_dict, dt=self.time_step*ms)
        
        regular_source = np.asarray(regular_source).astype(int)
        regular_target = np.asarray(regular_target).astype(int)
        
        self.regularsynapse.connect(i=regular_source, j=regular_target)
        
        regular_type = np.asarray(regular_type)
        
        self.regularsynapse.AMPA = 2 - regular_type
        self.regularsynapse.GABA = regular_type - 1
        self.regularsynapse.NMDA = 2 - regular_type
        
        self.regularsynapse.failure_rate = pfailinp
        self.regularsynapse.g_max = regular_gmax *nS
        
        self.regularsynapse.Gsyn_AMPA = self.Gsyn_AMPA
        self.regularsynapse.Gsyn_GABA = self.Gsyn_GABA
        self.regularsynapse.Gsyn_NMDA = self.Gsyn_NMDA
        
        self.regularspikemonitor = SpikeMonitor(self.regularinput)
        
        self.net.add(self.regularinput, self.regularsynapse, self.regularspikemonitor)
                   
    def run(self, t):
        
        if len(self.ta) > 0:
            ta = TimedArray(self.ta, dt=self.time_step*ms)
   
        self.net.run(t, report='text', report_period=10*second)
    
    def plot_V_t(self, ind, sim_dir, tlim=False, Vlim=False, V_Tlabel=False):
        
        if not os.path.isdir('{}/graph_V_t'.format(sim_dir)):
            os.mkdir('{}/graph_V_t'.format(sim_dir))
 
        
            
        curr_idc = 0
        NPfile = '{}/graph_V_t/V_t_{}_{}.png'.format(sim_dir, ind, curr_idc)
        
        while os.path.isfile(NPfile):
            curr_idc += 1
            NPfile = '{}/graph_V_t/V_t_{}_{}.png'.format(sim_dir, ind, curr_idc)
            
        fig, ax = subplots()
        
        V_mon = self.neuronmonitors['V']
        
        pos = np.where(V_mon.record==ind)[0][0]
        
        ax.plot(V_mon.t/ms, V_mon.V[pos]/mV, label='V')
        tspikes = self.spikemonitor.t[self.spikemonitor.i==pos]/ms
        
        for t in tspikes:
            ax.vlines(t, self.NeuPar[4, pos], 20)
        ax.set_xlabel('t (ms)')
        ax.set_ylabel('V (mV)')
        
        if type(tlim) is not bool:
            ax.set_xlim(tlim[0], tlim[1])
            
        if type(Vlim) is not bool:
            ax.set_ylim(Vlim[0], Vlim[1])
        if V_Tlabel:
            V_T = self.NeuPar[9, pos]
            if type(tlim) is not bool:
                t0, t1 = tlim
            else:
                t0 = 0
                t1 = V_mon.t[-1]/ms
            ax.hlines(V_T, t0, t1, label='V_T', linestyle='--', color='black')
        
        fig.legend()
          
        
        fig.savefig(NPfile)
       
    def plot_w_t(self, ind, sim_dir, tlim=False, wlim=False):
        
        if not os.path.isdir('{}/graph_w_t'.format(sim_dir)):
            os.mkdir('{}/graph_w_t'.format(sim_dir))
 
        
            
        curr_idc = 0
        NPfile = '{}/graph_w_t/w_t_{}_{}.png'.format(sim_dir, ind, curr_idc)
        
        while os.path.isfile(NPfile):
            curr_idc += 1
            NPfile = '{}/graph_w_t/w_t_{}_{}.png'.format(sim_dir, ind, curr_idc)
            
        fig, ax = subplots()
        
        w_mon = self.neuronmonitors['w']
        
        pos = np.where(w_mon.record==ind)[0][0]
        
        ax.plot(w_mon.t/ms, w_mon.w[pos]/pA)
        tspikes = self.spikemonitor.t[self.spikemonitor.i==pos]/ms
        
        for t in tspikes:
            ax.vlines(t, min(w_mon.w[pos]/pA), max(w_mon.w[pos]/pA), linestyle='--', color='black')
        ax.set_xlabel('t (ms)')
        ax.set_ylabel('w (pA)')
        
        if type(tlim) is not bool:
            ax.set_xlim(tlim[0], tlim[1])
            
        if type(wlim) is not bool:
            ax.set_ylim(wlim[0], wlim[1])
          
        
        fig.savefig(NPfile)



    def group_name_to_neuron_idc(self, name):
        
        idcs = []
        
        for analysand in name:
            grname, stripe = analysand
            
            if grname == 'all' and stripe == 'all':
            
                for gr in range(len(self.group_distr)):
                    for strp in range(len(self.group_distr[gr])):
                        idcs.extend(self.group_distr[gr][strp])

            else:
                if type(grname) is str:
                    group_idc = self.group_name_to_group_idc(grname)
                elif type(grname) is int:
                    group_idc = [grname,]
                elif type(grname) is list:
                    group_idc = grname
                    
                for group in group_idc:
                    idcs.extend(self.group_distr[group][stripe])

        return np.unique(np.asarray(idcs)).astype(int)

    def get_source_from_target(self, ind):
        
        _ = np.isin(self.target_arr, ind)
        
        return self.source_arr[_]

    def get_target_from_source(self, ind):
    
        _ = np.isin(self.source_arr, ind)
        
        return self.target_arr[_]
    
    def get_neurons_from_synapses(self, ind):
        
        return self.target_arr[ind], self.source_arr[ind]
        
    def group_name_to_synapses_idc(self, names, neuron_index=False):
        
        both = []
        
        for name in names:
            Nsyn = len(self.source_arr)
            
            if neuron_index:
                target_idc, source_idc = name
                
                if type(target_idc) is str and target_idc=='all':
                    target_idc = np.arange(Nsyn)
                if type(source_idc) is str and source_idc == 'all':
                    source_idc = np.arange(Nsyn)
                
            else:
                      
                target_name, target_stripe, source_name, source_stripe = name
    
                    
                target_idc = self.group_name_to_neuron_idc([[target_name, target_stripe],])
                source_idc = self.group_name_to_neuron_idc([[source_name, source_stripe],])
            
            Nsyn = len(self.source_arr)
            
            target_isin = np.arange(Nsyn)[np.isin(self.target_arr, target_idc)]
            source_isin = np.arange(Nsyn)[np.isin(self.source_arr, source_idc)]
            both.extend(np.intersect1d(target_isin, source_isin))
        
        both.sort()
        
        return np.asarray(both)

    def group_name_to_AMPAsynapses_idc(self, names, neuron_index=False):
        
        synidc = self.group_name_to_synapses_idc(names, neuron_index=neuron_index)
        synparidc = self.SynPar[0, synidc]
        AMPA = np.where(synparidc==1)[0]
           
        return synidc[AMPA]
    
    def group_name_to_GABAsynapses_idc(self, names, neuron_index=False):
        
        synidc = self.group_name_to_synapses_idc(names, neuron_index=neuron_index)
        synparidc = self.SynPar[0, synidc]
        GABA = np.where(synparidc==2)[0]
           
        return synidc[GABA]

    def group_name_to_NMDAsynapses_idc(self, names, neuron_index=False):
        
        synidc = self.group_name_to_synapses_idc(names, neuron_index=neuron_index)
        synparidc = self.SynPar[0, synidc]
        NMDA = np.where(synparidc==3)[0]
           
        return synidc[NMDA]
    
    def find_pCon(self, names, neuron_index=False):
        
        inhsyn = len(self.group_name_to_GABAsynapses_idc(names, neuron_index=neuron_index))
        excsyn = len(self.group_name_to_AMPAsynapses_idc(names, neuron_index=neuron_index))
        
        for name in names:
            if neuron_index:
                target_idc, source_idc = name
                
                if type(target_idc) is str and target_idc=='all':
                    Ntarget = self.NeuPar.shape[1]
                else:
                    Ntarget = len(target_idc)
                if type(source_idc) is str and source_idc == 'all':
                    Nsource = self.NeuPar.shape[1]
                else:
                    Nsource=len(source_idc)
                
            else:
                target_name, target_stripe, source_name, source_stripe = name
    
                    
                Ntarget = len(self.group_name_to_neuron_idc([[target_name, target_stripe],]))
                Nsource = len(self.group_name_to_neuron_idc([[source_name, source_stripe],]))
            
            totalsyn = Nsource*Ntarget
            
        return (inhsyn + excsyn)/totalsyn
    
    def find_inh_pCon(self, names, neuron_index=False):
        
        inhsyn = len(self.group_name_to_GABAsynapses_idc(names, neuron_index=neuron_index))
        
        totalsyn = 0
        for name in names:
            if neuron_index:
                target_idc, source_idc = name
                
                if type(target_idc) is str and target_idc=='all':
                    target_idc = np.arange(self.NeuPar.shape[1])
                if type(source_idc) is str and source_idc == 'all':
                    source_idc = np.arange(self.NeuPar.shape[1])
                
            else:
                target_name, target_stripe, source_name, source_stripe = name         
                target_idc = self.group_name_to_neuron_idc([[target_name, target_stripe],])
                source_idc = self.group_name_to_neuron_idc([[source_name, source_stripe],])
            
            inhidc=[]
            for col in range(len(self.group_distr[0])):
                inhidc.extend(list(self.group_name_to_neuron_idc([['IN', col]])))
            
            source_idc = source_idc[np.isin(source_idc, np.asarray(inhidc))]
            
            totalsyn += len(target_idc)*len(source_idc)
            
        if totalsyn == 0:
            inh_pCon = 'No inhibitory input'
        else:
            inh_pCon = inhsyn/totalsyn
        
        return inh_pCon

    def find_exc_pCon(self, names, neuron_index=False):
        
        excsyn = len(self.group_name_to_AMPAsynapses_idc(names, neuron_index=neuron_index))
        
        totalsyn = 0
        for name in names:
            if neuron_index:
                target_idc, source_idc = name
                
                if type(target_idc) is str and target_idc=='all':
                    target_idc = np.arange(self.NeuPar.shape[1])
                if type(source_idc) is str and source_idc == 'all':
                    source_idc = np.arange(self.NeuPar.shape[1])
                
            else:
                target_name, target_stripe, source_name, source_stripe = name         
                target_idc = self.group_name_to_neuron_idc([[target_name, target_stripe],])
                source_idc = self.group_name_to_neuron_idc([[source_name, source_stripe],])
            
            excidc=[]
            for col in range(len(self.group_distr[0])):
                excidc.extend(list(self.group_name_to_neuron_idc([['PC', col]])))
            
            source_idc = source_idc[np.isin(source_idc, np.asarray(excidc))]
            
            totalsyn += len(target_idc)*len(source_idc)
            
        if totalsyn == 0:
            exc_pCon = 'No excitatory input'
        else:
            exc_pCon = excsyn/totalsyn
        
        return exc_pCon

    def neuron_params(self, NeuronsSetup, sim_dir, bins=False):
        
        if not os.path.isdir('{}/NeuronParams'.format(sim_dir)):
            os.mkdir('{}/NeuronParams'.format(sim_dir))
 
        for neuron_idc in NeuronsSetup: 
            
            curr_idc = 0
            NPfile = '{}/NeuronParams/NeuronParams_{}.txt'.format(sim_dir, curr_idc)

            while os.path.isfile(NPfile):
                curr_idc += 1
                NPfile = '{}/NeuronParams/NeuronParams_{}.txt'.format(sim_dir, curr_idc)  
            
            
            with open(NPfile, 'a') as f:
                
                C_mean =  round(np.mean(self.NeuPar[0, neuron_idc]),3)
                C_std = round(np.std(self.NeuPar[0, neuron_idc]),3)
                C_min = round(np.min(self.NeuPar[0, neuron_idc]),3)
                C_max = round(np.max(self.NeuPar[0, neuron_idc]),3)
                
                g_L_mean =  round(np.mean(self.NeuPar[1, neuron_idc]),3)
                g_L_std = round(np.std(self.NeuPar[1, neuron_idc]),3)
                g_L_min = round(np.min(self.NeuPar[1, neuron_idc]),3)
                g_L_max = round(np.max(self.NeuPar[1, neuron_idc]),3)
                
                E_L_mean = round(np.mean(self.NeuPar[2, neuron_idc]),3)
                E_L_std = round(np.std(self.NeuPar[2, neuron_idc]),3)
                E_L_min = round(np.min(self.NeuPar[2, neuron_idc]),3)
                E_L_max = round(np.max(self.NeuPar[2, neuron_idc]),3)
                
                delta_T_mean =  round(np.mean(self.NeuPar[3, neuron_idc]),3)
                delta_T_std = round(np.std(self.NeuPar[3, neuron_idc]),3)
                delta_T_min = round(np.min(self.NeuPar[3, neuron_idc]),3)
                delta_T_max = round(np.max(self.NeuPar[3, neuron_idc]),3)
                
                V_up_mean =  round(np.mean(self.NeuPar[4, neuron_idc]),3)
                V_up_std = round(np.std(self.NeuPar[4, neuron_idc]),3)
                V_up_min = round(np.min(self.NeuPar[4, neuron_idc]),3)
                V_up_max = round(np.max(self.NeuPar[4, neuron_idc]),3)
                               
                tau_w_mean =  round(np.mean(self.NeuPar[5, neuron_idc]),3)
                tau_w_std = round(np.std(self.NeuPar[5, neuron_idc]),3)
                tau_w_min = round(np.min(self.NeuPar[5, neuron_idc]),3)
                tau_w_max = round(np.max(self.NeuPar[5, neuron_idc]),3)
                           
                
                b_mean =  round(np.mean(self.NeuPar[7, neuron_idc]),3)
                b_std = round(np.std(self.NeuPar[7, neuron_idc]),3)
                b_min = round(np.min(self.NeuPar[7, neuron_idc]),3)
                b_max = round(np.max(self.NeuPar[7, neuron_idc]),3)
                
                V_r_mean =  round(np.mean(self.NeuPar[8, neuron_idc]),3)
                V_r_std = round(np.std(self.NeuPar[8, neuron_idc]),3)
                V_r_min = round(np.min(self.NeuPar[8, neuron_idc]),3)
                V_r_max = round(np.max(self.NeuPar[8, neuron_idc]),3)
                
                V_T_mean =  round(np.mean(self.NeuPar[9, neuron_idc]),3)
                V_T_std = round(np.std(self.NeuPar[9, neuron_idc]),3)
                V_T_min = round(np.min(self.NeuPar[9, neuron_idc]),3)
                V_T_max = round(np.max(self.NeuPar[9, neuron_idc]),3)
                
                I_ref_mean =  round(np.mean(self.NeuPar[10, neuron_idc]),3)
                I_ref_std = round(np.std(self.NeuPar[10, neuron_idc]),3)
                I_ref_min = round(np.min(self.NeuPar[10, neuron_idc]),3)
                I_ref_max = round(np.max(self.NeuPar[10, neuron_idc]),3)
                
                V_dep_mean =  round(np.mean(self.NeuPar[11, neuron_idc]),3)
                V_dep_std = round(np.std(self.NeuPar[11, neuron_idc]),3)
                V_dep_min = round(np.min(self.NeuPar[11, neuron_idc]),3)
                V_dep_max = round(np.max(self.NeuPar[11, neuron_idc]),3)
                
                tau_m = self.NeuPar[0, neuron_idc]/self.NeuPar[1, neuron_idc]
                tau_m_mean =  round(np.mean(tau_m),3)
                tau_m_std = round(np.std(tau_m),3)
                tau_m_min = round(np.min(tau_m),3)
                tau_m_max = round(np.max(tau_m),3)
                
                I_SN = self.I_SN[neuron_idc]/pA
                I_SN_mean = round(np.mean(I_SN), 3)
                I_SN_std = round(np.std(I_SN),3)
                I_SN_min = round(np.min(I_SN),3)
                I_SN_max = round(np.max(I_SN),3)
                
                text = "Neurons:{}\n\nmean += std (min, max)\n\n".format(neuron_idc)
                text+="C = {}+-{} (min: {}, max:  {}) pF\n".format(C_mean, C_std, C_min, C_max)
                text+="g_L = {} +-{} (min: {}, max:  {}) nS\n".format(g_L_mean, g_L_std, g_L_min, g_L_max)
                text+="E_L = {} +-{} (min: {}, max:  {}) mV\n\n".format(E_L_mean, E_L_std, E_L_min, E_L_max)
                
                text+="delta_T = {} +-{} (min: {}, max:  {}) mV\n".format(delta_T_mean, delta_T_std, delta_T_min, delta_T_max)
                text+="V_T = {} +-{} (min: {}, max:  {}) mV\n".format(V_T_mean, V_T_std, V_T_min, V_T_max)
                text+="V_up = {} +-{} (min: {}, max:  {}) mV\n".format(V_up_mean, V_up_std, V_up_min, V_up_max)
                text+="V_r = {} +-{} (min: {}, max:  {}) mV\n\n".format(V_r_mean, V_r_std, V_r_min, V_r_max)
                
                text+="b = {} +-{} (min: {}, max:  {}) pA\n".format(b_mean, b_std, b_min, b_max)
                text+="tau_w = {} +-{} (min: {}, max:  {}) ms\n".format(tau_w_mean, tau_w_std, tau_w_min, tau_w_max)
                text+="tau_m = {} +-{} (min: {}, max:  {}) ms\n".format(tau_m_mean, tau_m_std, tau_m_min, tau_m_max)
    
                text+="I_SN = {} +-{} (min: {}, max:  {}) pA\n".format(I_SN_mean, I_SN_std, I_SN_min,I_SN_max)
                text+="V_SN = V_T = {} +-{} (min: {}, max:  {}) mV\n\n".format(V_T_mean, V_T_std, V_T_min, V_T_max)
                
                text+="V_dep = V_r = {} +-{} (min: {}, max:  {}) mV\n".format(V_dep_mean, V_dep_std, V_dep_min, V_dep_max)
                text+="I_ref = {} +-{} (min: {}, max:  {}) pA\n\n".format(I_ref_mean, I_ref_std, I_ref_min, I_ref_max)
                text+="-"*50 + '\n'
    
                
                print(text, file=f)
                
             
                
                    
                if not bins:
                    bins = int(round(np.sqrt(len(neuron_idc)), 0))
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[0, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('C (pF)')
                fig.savefig("{}/NeuronParams/C_hist_{}.png".format(sim_dir, curr_idc))
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[0, neuron_idc])
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[0, neuron_idc], bins=bins_lims, density=True)
                ax.set_ylabel('Density')
                ax.set_xlabel('C (pF)')
                fig.savefig("{}/NeuronParams/C_histdens_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('C (pF)')
                fig.savefig("{}/NeuronParams/C_hist_plot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('C (pF)')
                ax.set_xlim(np.mean(self.NeuPar[0, neuron_idc]) - 1*np.std(self.NeuPar[0, neuron_idc]), np.mean(self.NeuPar[0, neuron_idc]) + 1*np.std(self.NeuPar[0, neuron_idc]))
                fig.savefig("{}/NeuronParams/C_hist_plot1_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('C (pF)')
                ax.set_xlim(np.mean(self.NeuPar[0, neuron_idc]) - 2*np.std(self.NeuPar[0, neuron_idc]), np.mean(self.NeuPar[0, neuron_idc]) + 2*np.std(self.NeuPar[0, neuron_idc]))
                fig.savefig("{}/NeuronParams/C_hist_plot2_{}.png".format(sim_dir, curr_idc))
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[1, neuron_idc])
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[1, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('g_L (nS)')
                fig.savefig("{}/NeuronParams/g_L_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[1, neuron_idc], bins=bins_lims, density=True)
                ax.set_ylabel('Density')
                ax.set_xlabel('g_L (nS)')
                fig.savefig("{}/NeuronParams/g_L_histdens_{}.png".format(sim_dir, curr_idc))
                
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('g_L (nS)')
                fig.savefig("{}/NeuronParams/g_L_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('g_L (nS)')
                ax.set_xlim(0, np.mean(self.NeuPar[1, neuron_idc]))
                fig.savefig("{}/NeuronParams/g_L_histplot0_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('g_L (nS)')
                ax.set_xlim(0, np.mean(self.NeuPar[1, neuron_idc]) + 1 * np.std(self.NeuPar[1, neuron_idc]))
                fig.savefig("{}/NeuronParams/g_L_histplot1_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('g_L (nS)')
                ax.set_xlim(0, np.mean(self.NeuPar[1, neuron_idc]) + 2 * np.std(self.NeuPar[1, neuron_idc]))
                fig.savefig("{}/NeuronParams/g_L_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[2, neuron_idc])
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[2, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('E_L (mV)')
                fig.savefig("{}/NeuronParams/E_L_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[2, neuron_idc], bins=bins_lims, density=True)
                ax.set_ylabel('Density')
                ax.set_xlabel('E_L (mV)')
                fig.savefig("{}/NeuronParams/E_L_histdens_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('E_L (mV)')
                fig.savefig("{}/NeuronParams/E_L_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('E_L (mV)')
                ax.set_xlim(np.mean(self.NeuPar[2, neuron_idc]) - 1*np.std(self.NeuPar[2, neuron_idc]), np.mean(self.NeuPar[2, neuron_idc]) + 1*np.std(self.NeuPar[2, neuron_idc]))
                fig.savefig("{}/NeuronParams/E_L_histplot1_{}.png".format(sim_dir, curr_idc))
                
                ig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('E_L (mV)')
                ax.set_xlim(np.mean(self.NeuPar[2, neuron_idc]) - 2*np.std(self.NeuPar[2, neuron_idc]), np.mean(self.NeuPar[2, neuron_idc]) + 2*np.std(self.NeuPar[2, neuron_idc]))
                fig.savefig("{}/NeuronParams/E_L_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[3, neuron_idc])
        
                fig, ax = subplots()
                ax.hist(self.NeuPar[3, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('delta_T (mV)')
                fig.savefig("{}/NeuronParams/delta_T_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[3, neuron_idc], bins=bins_lims, density=True)
                ax.set_ylabel('Density')
                ax.set_xlabel('delta_T (mV)')
                fig.savefig("{}/NeuronParams/delta_T_histdens_{}.png".format(sim_dir, curr_idc))
                
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('delta_T (mV)')
                fig.savefig("{}/NeuronParams/delta_T_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('delta_T (mV)')
                ax.set_xlim(0, np.mean(self.NeuPar[3, neuron_idc]))
                fig.savefig("{}/NeuronParams/delta_T_histplot0_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('delta_T (mV)')
                ax.set_xlim(0, np.mean(self.NeuPar[3, neuron_idc]) + 1 * np.std(self.NeuPar[3, neuron_idc]))
                fig.savefig("{}/NeuronParams/delta_T_histplot1_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('delta_T (mV)')
                ax.set_xlim(0, np.mean(self.NeuPar[3, neuron_idc]) + 2 * np.std(self.NeuPar[3, neuron_idc]))
                fig.savefig("{}/NeuronParams/delta_T_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[4, neuron_idc])
                                
                fig, ax = subplots()
                ax.hist(self.NeuPar[4, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                fig.savefig("{}/NeuronParams/V_up_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[4, neuron_idc], bins=bins_lims, density=True)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_up (mV)')
                fig.savefig("{}/NeuronParams/V_up_histdens_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_up (mV)')
                fig.savefig("{}/NeuronParams/V_up_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_up (mV)')
                ax.set_xlim(np.mean(self.NeuPar[4, neuron_idc]) - 1*np.std(self.NeuPar[4, neuron_idc]), np.mean(self.NeuPar[4, neuron_idc]) + 1*np.std(self.NeuPar[4, neuron_idc]))
                fig.savefig("{}/NeuronParams/V_up_histplot1_{}.png".format(sim_dir, curr_idc))
                
                ig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_up (mV)')
                ax.set_xlim(np.mean(self.NeuPar[4, neuron_idc]) - 2*np.std(self.NeuPar[4, neuron_idc]), np.mean(self.NeuPar[4, neuron_idc]) + 2*np.std(self.NeuPar[4, neuron_idc]))
                fig.savefig("{}/NeuronParams/V_up_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[5, neuron_idc])
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[5, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('tau_w (ms)')
                fig.savefig("{}/NeuronParams/tau_w_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('tau_w (ms)')
                fig.savefig("{}/NeuronParams/tau_w_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('tau_w (ms)')
                ax.set_xlim(0, np.mean(self.NeuPar[5, neuron_idc]))
                fig.savefig("{}/NeuronParams/tau_w_histplot0_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('tau_w (ms)')
                ax.set_xlim(0, np.mean(self.NeuPar[5, neuron_idc]) + 1 * np.std(self.NeuPar[5, neuron_idc]))
                fig.savefig("{}/NeuronParams/tau_w_histplot1_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('tau_w (ms)')
                ax.set_xlim(0, np.mean(self.NeuPar[5, neuron_idc]) + 2 * np.std(self.NeuPar[5, neuron_idc]))
                fig.savefig("{}/NeuronParams/tau_w_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[7, neuron_idc])
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[7, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('b (pA)')
                fig.savefig("{}/NeuronParams/b_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('b (pA)')
                fig.savefig("{}/NeuronParams/b_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('b (pA)')
                ax.set_xlim(0, np.mean(self.NeuPar[7, neuron_idc]))
                fig.savefig("{}/NeuronParams/b_histplot0_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('b (pA)')
                ax.set_xlim(0, np.mean(self.NeuPar[7, neuron_idc]) + 1 * np.std(self.NeuPar[7, neuron_idc]))
                fig.savefig("{}/NeuronParams/b_histplot1_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('b (pA)')
                ax.set_xlim(0, np.mean(self.NeuPar[7, neuron_idc]) + 2 * np.std(self.NeuPar[7, neuron_idc]))
                fig.savefig("{}/NeuronParams/b_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[8, neuron_idc])
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[8, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('V_r (mV)')
                fig.savefig("{}/NeuronParams/V_r_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[8, neuron_idc], bins=bins_lims, density=True)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_r (mV)')
                fig.savefig("{}/NeuronParams/V_r_histdens_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_r (mV)')
                fig.savefig("{}/NeuronParams/V_r_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_r (mV)')
                ax.set_xlim(np.mean(self.NeuPar[8, neuron_idc]) - 1*np.std(self.NeuPar[8, neuron_idc]), np.mean(self.NeuPar[8, neuron_idc]) + 1*np.std(self.NeuPar[8, neuron_idc]))
                fig.savefig("{}/NeuronParams/V_r_histplot1_{}.png".format(sim_dir, curr_idc))
                
                ig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_r (mV)')
                ax.set_xlim(np.mean(self.NeuPar[8, neuron_idc]) - 2*np.std(self.NeuPar[8, neuron_idc]), np.mean(self.NeuPar[8, neuron_idc]) + 2*np.std(self.NeuPar[8, neuron_idc]))
                fig.savefig("{}/NeuronParams/V_r_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[9, neuron_idc])
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[9, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('V_T (mV)')
                fig.savefig("{}/NeuronParams/V_T_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[9, neuron_idc], bins=bins_lims, density=True)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_T (mV)')
                fig.savefig("{}/NeuronParams/V_T_histdens_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_T (mV)')
                fig.savefig("{}/NeuronParams/V_T_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_T (mV)')
                ax.set_xlim(np.mean(self.NeuPar[9, neuron_idc]) - 1*np.std(self.NeuPar[9, neuron_idc]), np.mean(self.NeuPar[9, neuron_idc]) + 1*np.std(self.NeuPar[9, neuron_idc]))
                fig.savefig("{}/NeuronParams/V_T_histplot1_{}.png".format(sim_dir, curr_idc))
                
                ig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('V_T (mV)')
                ax.set_xlim(np.mean(self.NeuPar[9, neuron_idc]) - 2*np.std(self.NeuPar[9, neuron_idc]), np.mean(self.NeuPar[9, neuron_idc]) + 2*np.std(self.NeuPar[9, neuron_idc]))
                fig.savefig("{}/NeuronParams/V_T_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                
                
                
                bins_lims, dens_lims, dens_arr = build_hist(self.NeuPar[10, neuron_idc])
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[10, neuron_idc], bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('I_ref (pA)')
                fig.savefig("{}/NeuronParams/I_ref_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(self.NeuPar[10, neuron_idc], bins=bins_lims, density=True)
                ax.set_ylabel('Density')
                ax.set_xlabel('I_ref (pA)')
                fig.savefig("{}/NeuronParams/I_ref_histdens_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('I_ref (pA)')
                fig.savefig("{}/NeuronParams/I_ref_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('I_ref (pA)')
                ax.set_xlim(np.mean(self.NeuPar[10, neuron_idc]) - 1*np.std(self.NeuPar[10, neuron_idc]), np.mean(self.NeuPar[10, neuron_idc]) + 1*np.std(self.NeuPar[10, neuron_idc]))
                fig.savefig("{}/NeuronParams/I_ref_histplot1_{}.png".format(sim_dir, curr_idc))
                
                ig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('I_ref (pA)')
                ax.set_xlim(np.mean(self.NeuPar[10, neuron_idc]) - 2*np.std(self.NeuPar[10, neuron_idc]), np.mean(self.NeuPar[10, neuron_idc]) + 2*np.std(self.NeuPar[10, neuron_idc]))
                fig.savefig("{}/NeuronParams/I_ref_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                
                
                bins_lims, dens_lims, dens_arr = build_hist(tau_m)
                                
                fig, ax = subplots()
                ax.hist(tau_m, bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('tau_m (ms)')
                fig.savefig("{}/NeuronParams/tau_m_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('tau_m (ms)')
                fig.savefig("{}/NeuronParams/tau_m_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('tau_m (ms)')
                ax.set_xlim(0, np.mean(tau_m))
                fig.savefig("{}/NeuronParams/tau_m_histplot0_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('tau_m (ms)')
                ax.set_xlim(0, np.mean(tau_m) + 1 * np.std(tau_m))
                fig.savefig("{}/NeuronParams/tau_m_histplot1_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('tau_m (ms)')
                ax.set_xlim(0, np.mean(tau_m) + 2 * np.std(tau_m))
                fig.savefig("{}/NeuronParams/tau_m_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                bins_lims, dens_lims, dens_arr = build_hist(I_SN)
                
                fig, ax = subplots()
                ax.hist(I_SN, bins=bins)
                ax.set_ylabel('Number of neurons')
                ax.set_xlabel('I_SN (pA)')
                fig.savefig("{}/NeuronParams/I_SN_hist_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(tau_m, bins=bins_lims, density=True)
                ax.set_ylabel('Density')
                ax.set_xlabel('I_SN (pA)')
                fig.savefig("{}/NeuronParams/I_SN_histdens_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('I_SN (pA)')
                fig.savefig("{}/NeuronParams/I_SN_histplot_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('I_SN (pA)')
                ax.set_xlim(np.mean(I_SN) - 1*np.std(I_SN), np.mean(I_SN) + 1*np.std(I_SN))
                fig.savefig("{}/NeuronParams/I_SN_histplot1_{}.png".format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_ylabel('Density')
                ax.set_xlabel('I_SN (pA)')
                ax.set_xlim(np.mean(I_SN) - 2*np.std(I_SN), np.mean(I_SN) + 2*np.std(I_SN))
                fig.savefig("{}/NeuronParams/I_SN_histplot2_{}.png".format(sim_dir, curr_idc))
                
                
                
                
                
                
                for ind in neuron_idc:
                    group, stripe = self.neuron_idc_to_group_name(ind)
                    text = "Neuron {}, Group: {}, Stripe: {}\n".format(ind, group, stripe)
        
                    text+= "C = {} pF\n".format(round(self.NeuPar[0, ind], 3))
                    text+= "g_L = {} nS\n".format(round(self.NeuPar[1, ind] ,3))
                    text+= "E_L = {} mV\n\n".format(round(self.NeuPar[2, ind] ,3))
                    
                    text+= "delta_T = {} mV\n".format(round(self.NeuPar[3, ind] ,3))
                    text+= "V_T = {} mV\n".format(round(self.NeuPar[9, ind] ,3))
                    text+= "V_up = {} mV\n".format(round(self.NeuPar[4, ind] ,3))
                    text+= "V_r = {} mV\n\n".format(round(self.NeuPar[8, ind] ,3))
                    
                    text+= "b = {} pA\n".format(round(self.NeuPar[7, ind] ,3))
                    text+= "tau_w = {} ms\n".format(round(self.NeuPar[5, ind] ,3))
                    text+= "tau_m = {} ms\n\n".format(round(self.NeuPar[0, ind]/self.NeuPar[1, ind] ,3))
        
                    text+= "I_SN = {} pA\n".format(round(self.I_SN[ind]/pA ,3))
                    text+= "V_SN = V_T = {} mV\n\n".format(round(self.NeuPar[9, ind] ,3))
                    
                    text+= "V_dep = V_r = {} mV\n".format(round(self.NeuPar[11, ind] ,3))
                    text+= "I_ref = {} pA\n\n".format(round(self.NeuPar[10, ind] ,3))
                    
                    text+="-"*50 + "\n"
        
        
                    print(text, file=f)
                
    def groups_params(self, GroupSetup, sim_dir, bins=False, histograms=True):
        
        if not os.path.isdir('{}/GroupParams'.format(sim_dir)):
            os.mkdir('{}/GroupParams'.format(sim_dir))
        
            
        for groups in GroupSetup: 

            curr_idc = 0
            NPfile = '{}/GroupParams/GroupParams_{}.txt'.format(sim_dir, curr_idc)
            
            while os.path.isfile(NPfile):
                curr_idc += 1
                NPfile = '{}/GroupParams/GroupParams_{}.txt'.format(sim_dir, curr_idc)
             
            neuron_idc = np.asarray(self.group_name_to_neuron_idc(groups))
            
            group_stripe = []
           
            for group in groups:
                grname, stripe = group
                group_stripe.append('{} (stripe {})'.format(grname, stripe))
            
            if len(neuron_idc)>0:
                groups_stripes = ', '.join(group_stripe)
                
                C_mean = round(np.mean(self.NeuPar[0, neuron_idc]),3)
                C_std = round(np.std(self.NeuPar[0, neuron_idc]),3)
                C_min = round(np.min(self.NeuPar[0, neuron_idc]),3)
                C_max = round(np.max(self.NeuPar[0, neuron_idc]),3)
                
                g_L_mean =  round(np.mean(self.NeuPar[1, neuron_idc]),3)
                g_L_std = round(np.std(self.NeuPar[1, neuron_idc]),3)
                g_L_min = round(np.min(self.NeuPar[1, neuron_idc]),3)
                g_L_max = round(np.max(self.NeuPar[1, neuron_idc]),3)
                
                E_L_mean = round(np.mean(self.NeuPar[2, neuron_idc]),3)
                E_L_std = round(np.std(self.NeuPar[2, neuron_idc]),3)
                E_L_min = round(np.min(self.NeuPar[2, neuron_idc]),3)
                E_L_max = round(np.max(self.NeuPar[2, neuron_idc]),3)
                
                delta_T_mean = round(np.mean(self.NeuPar[3, neuron_idc]),3)
                delta_T_std = round(np.std(self.NeuPar[3, neuron_idc]),3)
                delta_T_min = round(np.min(self.NeuPar[3, neuron_idc]),3)
                delta_T_max = round(np.max(self.NeuPar[3, neuron_idc]),3)
                
                V_up_mean = round(np.mean(self.NeuPar[4, neuron_idc]),3)
                V_up_std = round(np.std(self.NeuPar[4, neuron_idc]),3)
                V_up_min = round(np.min(self.NeuPar[4, neuron_idc]),3)
                V_up_max = round(np.max(self.NeuPar[4, neuron_idc]),3)
                               
                tau_w_mean = round(np.mean(self.NeuPar[5, neuron_idc]),3)
                tau_w_std = round(np.std(self.NeuPar[5, neuron_idc]),3)
                tau_w_min = round(np.min(self.NeuPar[5, neuron_idc]),3)
                tau_w_max = round(np.max(self.NeuPar[5, neuron_idc]),3)
                           
                
                b_mean = round(np.mean(self.NeuPar[7, neuron_idc]),3)
                b_std = round(np.std(self.NeuPar[7, neuron_idc]),3)
                b_min = round(np.min(self.NeuPar[7, neuron_idc]),3)
                b_max = round(np.max(self.NeuPar[7, neuron_idc]),3)
                
                V_r_mean = round(np.mean(self.NeuPar[8, neuron_idc]),3)
                V_r_std = round(np.std(self.NeuPar[8, neuron_idc]),3)
                V_r_min = round(np.min(self.NeuPar[8, neuron_idc]),3)
                V_r_max = round(np.max(self.NeuPar[8, neuron_idc]),3)
                
                V_T_mean = round(np.mean(self.NeuPar[9, neuron_idc]),3)
                V_T_std = round(np.std(self.NeuPar[9, neuron_idc]),3)
                V_T_min = round(np.min(self.NeuPar[9, neuron_idc]),3)
                V_T_max = round(np.max(self.NeuPar[9, neuron_idc]),3)
                
                I_ref_mean = round(np.mean(self.NeuPar[10, neuron_idc]),3)
                I_ref_std = round(np.std(self.NeuPar[10, neuron_idc]),3)
                I_ref_min = round(np.min(self.NeuPar[10, neuron_idc]),3)
                I_ref_max = round(np.max(self.NeuPar[10, neuron_idc]),3)
                
                V_dep_mean = round(np.mean(self.NeuPar[11, neuron_idc]),3)
                V_dep_std = round(np.std(self.NeuPar[11, neuron_idc]),3)
                V_dep_min = round(np.min(self.NeuPar[11, neuron_idc]),3)
                V_dep_max = round(np.max(self.NeuPar[11, neuron_idc]),3)
                
                tau_m = self.NeuPar[0, neuron_idc]/self.NeuPar[1, neuron_idc]
                
                tau_m_mean = round(np.mean(tau_m),3)
                tau_m_std = round(np.std(tau_m),3)
                tau_m_min = round(np.min(tau_m),3)
                tau_m_max = round(np.max(tau_m),3)
                
                I_SN = self.I_SN[neuron_idc]/pA
                I_SN_mean = round(np.mean(I_SN), 3)
                I_SN_std = round(np.std(I_SN),3)
                I_SN_min = round(np.min(I_SN),3)
                I_SN_max = round(np.max(I_SN),3)
                
                if not bins:
                    bins = int(round(np.sqrt(len(neuron_idc)), 0))
                
                if histograms:
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[0, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('C (pF)')
                    fig.savefig("{}/GroupParams/C_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[1, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('g_L (nS)')
                    fig.savefig("{}/GroupParams/g_L_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[2, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('E_L (mV)')
                    fig.savefig("{}/GroupParams/E_L_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[3, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('delta_T (ms)')
                    fig.savefig("{}/GroupParams/delta_T_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[4, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    fig.savefig("{}/GroupParams/V_up_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[5, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('tau_w (ms)')
                    fig.savefig("{}/GroupParams/tau_w_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[7, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    fig.savefig("{}/GroupParams/b_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[8, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('V_r (mV)')
                    fig.savefig("{}/GroupParams/V_r_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[9, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('V_T (mV)')
                    fig.savefig("{}/GroupParams/V_T_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[10, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('I_ref (pA)')
                    fig.savefig("{}/GroupParams/I_ref_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(self.NeuPar[11, neuron_idc], bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('V_dep (mV)')
                    fig.savefig("{}/GroupParams/V_dep_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(tau_m, bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('tau_m (ms)')
                    fig.savefig("{}/GroupParams/tau_m_hist_{}.png".format(sim_dir, curr_idc))
                    
                    fig, ax = subplots()
                    ax.hist(I_SN, bins=bins)
                    ax.set_ylabel('Number of neurons')
                    ax.set_xlabel('I_SN (pA)')
                    fig.savefig("{}/GroupParams/I_SN_hist_{}.png".format(sim_dir, curr_idc))
                
                
                     
                with open(NPfile, 'a') as f:
                    
                    text = "Group:{}\n\nNumber of neurons : {}\n\nMean (+- std) (min, max)\n".format(groups_stripes, len(neuron_idc))
                    text+="C = {}+-{} (min: {}, max:  {}) pF\n".format(C_mean, C_std, C_min, C_max)
                    text+="g_L = {} +-{} (min: {}, max:  {}) nS\n".format(g_L_mean, g_L_std, g_L_min, g_L_max)
                    text+="E_L = {} +-{} (min: {}, max:  {}) mV\n\n".format(E_L_mean, E_L_std, E_L_min, E_L_max)
                    
                    text+="delta_T = {} +-{} (min: {}, max:  {}) mV\n".format(delta_T_mean, delta_T_std, delta_T_min, delta_T_max)
                    text+="V_T = {} +-{} (min: {}, max:  {}) mV\n".format(V_T_mean, V_T_std, V_T_min, V_T_max)
                    text+="V_up = {} +-{} (min: {}, max:  {}) mV\n".format(V_up_mean, V_up_std, V_up_min, V_up_max)
                    text+="V_r = {} +-{} (min: {}, max:  {}) mV\n\n".format(V_r_mean, V_r_std, V_r_min, V_r_max)
                    
                    text+="b = {} +-{} (min: {}, max:  {}) pA\n".format(b_mean, b_std, b_min, b_max)
                    text+="tau_w = {} +-{} (min: {}, max:  {}) ms\n".format(tau_w_mean, tau_w_std, tau_w_min, tau_w_max)
                    text+="tau_m = {} +-{} (min: {}, max:  {}) ms\n".format(tau_m_mean, tau_m_std, tau_m_min, tau_m_max)
        
                    text+="I_SN = {} +-{} (min: {}, max:  {}) pA\n".format(I_SN_mean, I_SN_std, I_SN_min,I_SN_max)
                    text+="V_SN = V_T = {} +-{} (min: {}, max:  {}) mV\n\n".format(V_T_mean, V_T_std, V_T_min, V_T_max)
                    
                    text+="V_dep = V_r = {} +-{} (min: {}, max:  {}) mV\n".format(V_dep_mean, V_dep_std, V_dep_min, V_dep_max)
                    text+="I_ref = {} +-{} (min: {}, max:  {}) pA\n\n".format(I_ref_mean, I_ref_std, I_ref_min, I_ref_max)
                    text+="-"*50 + '\n'
                    print(text, file=f)
            
            else:
                with open(NPfile, 'a') as f:
                    
                    text = """
                    Group(s): {}
                    
                    No neurons
        
                    -----------------------------------
                    """.format(groups_stripes)
                    
                    print(text, file=f)

    def synapses_params(self, GroupSetup, sim_dir, bins=False, neuron_index=False):

        if not os.path.isdir('{}/SynapsesParams'.format(sim_dir)):
            os.mkdir('{}/SynapsesParams'.format(sim_dir))
        
            
        for groups in GroupSetup:    
            
            curr_idc = 0
            NPfile = '{}/SynapsesParams/SynapsesParams_{}.txt'.format(sim_dir, curr_idc)
            while os.path.isfile(NPfile):
                curr_idc += 1
                NPfile = '{}/SynapsesParams/SynapsesParams_{}.txt'.format(sim_dir, curr_idc)
            
            
            syn_idc = np.asarray(self.group_name_to_synapses_idc(groups, neuron_index=neuron_index))    
 
            group_stripe = []
           
            Ntargetxsource = 0
            for group in groups:
                if neuron_index:
                    target_idc, source_idc = group
                    
                    if type(target_idc) is str and target_idc=='all':
                        Ntarget = self.NeuPar.shape[1]
                    else:
                        Ntarget = len(target_idc)
                    
                    if type(source_idc) is str and source_idc=='all':
                        Nsource = self.NeuPar.shape[1]
                    else:
                        Nsource = len(source_idc)
                    
                    group_stripe.append('From {}\n\nto {}'.format(source_idc, target_idc))
                    Ntargetxsource += Ntarget*Nsource
                
                else:
                    grname_target, stripe_target, grname_source, stripe_source = group
                    group_stripe.append('From {} (stripe {}) to {} (stripe {})'.format(grname_source, stripe_source, grname_target, stripe_target))
                    Ntarget = len(self.group_name_to_neuron_idc([[grname_target, stripe_target]]))
                    Nsource = len(self.group_name_to_neuron_idc([[grname_source, stripe_source]]))
                    Ntargetxsource += Ntarget*Nsource                
            
            groups_stripes = ', '.join(group_stripe)
            
            SynP = self.SynPar[:, syn_idc]
            stypes = SynP[0, :]
            
            AMPA_idc = np.where(stypes==1)[0]
            NMDA_idc = np.where(stypes==3)[0]
            GABA_idc = np.where(stypes==2)[0]

            
            with open(NPfile, 'a') as f:
                f.write('Group(s): {}\n\n'.format(groups_stripes))
                          
            if len(AMPA_idc)> 0:
                use = SynP[1, AMPA_idc]
                tau_rec = SynP[2, AMPA_idc]
                tau_fac = SynP[3, AMPA_idc]
                g_max = SynP[4, AMPA_idc]
                dtax = SynP[5, AMPA_idc]
                p_fail = SynP[6, AMPA_idc]          
    
                if not bins:
                    bins = int(round(np.sqrt(len(AMPA_idc)), 0))
                    
                fig, ax = subplots()
                ax.hist(use, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('U')
                fig.savefig('{}/SynapsesParams/AMPA_U_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(tau_rec, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('tau_rec (ms)')
                fig.savefig('{}/SynapsesParams/AMPA_tau_rec_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(tau_fac, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('tau_fac (ms)')
                fig.savefig('{}/SynapsesParams/AMPA_tau_fac_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(g_max, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('g_max (nS)')
                fig.savefig('{}/SynapsesParams/AMPA_g_max_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(dtax, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('tau_D (ms)')
                fig.savefig('{}/SynapsesParams/AMPA_tau_D_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(p_fail, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('p_fail')
                fig.savefig('{}/SynapsesParams/AMPA_p_fail_hist_{}.png'.format(sim_dir, curr_idc))
                            
                    
                
    
                with open(NPfile, 'a') as f:
                    
                    text = "AMPA synapses\n\n"
                    
                    text+="Number of synapses: {}\n".format(len(AMPA_idc))
                    text+="Mean pCon: {}\n".format(len(AMPA_idc)/(Ntargetxsource))
                    
                    text+="Mean +- std (min, max)\n\n"
        
                    text+="U = {} +- {} (min: {}, max: {})\n".format(mean(use), std(use), min(use), max(use))
                    text+="tau_rec = {} +- {} (min: {}, max: {}) ms\n".format(mean(tau_rec), std(tau_rec), min(tau_rec), max(tau_rec))
                    text+="tau_fac = {} +- {} (min: {}, max: {}) ms\n\n".format(mean(tau_fac), std(tau_fac), min(tau_fac), max(tau_fac))
                    
                    text+="g_max = {} +- {} (min: {}, max: {}) nS\n".format(mean(g_max), std(g_max), min(g_max), max(g_max))
                    text+="delay = {} +- {} (min: {}, max: {}) ms\n".format(mean(dtax), std(dtax), min(dtax), max(dtax))
                    text+="p_fail = {} +- {}(min: {}, max: {})\n".format(mean(p_fail), std(p_fail), min(p_fail), max(p_fail))
                    text+='-'*50 + '\n'
                    
                    print(text, file=f)
            
            if len(NMDA_idc)> 0:
                use = SynP[1, NMDA_idc]
                tau_rec = SynP[2, NMDA_idc]
                tau_fac = SynP[3, NMDA_idc]
                g_max = SynP[4, NMDA_idc]
                dtax = SynP[5, NMDA_idc]
                p_fail = SynP[6, NMDA_idc]
                
                if not bins:
                    bins = int(round(np.sqrt(len(NMDA_idc)), 0))
                    
                fig, ax = subplots()
                ax.hist(use, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('U')
                fig.savefig('{}/SynapsesParams/NMDA_U_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(tau_rec, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('tau_rec (ms)')
                fig.savefig('{}/SynapsesParams/NMDA_tau_rec_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(tau_fac, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('tau_fac (ms)')
                fig.savefig('{}/SynapsesParams/NMDA_tau_fac_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(g_max, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('g_max (nS)')
                fig.savefig('{}/SynapsesParams/NMDA_g_max_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(dtax, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('tau_D (ms)')
                fig.savefig('{}/SynapsesParams/NMDA_tau_D_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(p_fail, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('p_fail')
                fig.savefig('{}/SynapsesParams/NMDA_p_fail_hist_{}.png'.format(sim_dir, curr_idc))
                            
                    
                
    
            
                with open(NPfile, 'a') as f:
                    
                  text = "NMDA synapses\n\n"
                  
                  text+="Number of synapses: {}\n".format(len(NMDA_idc))
                  text+="Mean pCon: {}\n".format(len(NMDA_idc)/(Ntargetxsource))
                  
                  text+="Mean +- std (min, max)\n\n"
      
                  text+="U = {} +- {} (min: {}, max: {})\n".format(mean(use), std(use), min(use), max(use))
                  text+="tau_rec = {} +- {} (min: {}, max: {}) ms\n".format(mean(tau_rec), std(tau_rec), min(tau_rec), max(tau_rec))
                  text+="tau_fac = {} +- {} (min: {}, max: {}) ms\n\n".format(mean(tau_fac), std(tau_fac), min(tau_fac), max(tau_fac))
                  
                  text+="g_max = {} +- {} (min: {}, max: {}) nS\n".format(mean(g_max), std(g_max), min(g_max), max(g_max))
                  text+="delay = {} +- {} (min: {}, max: {}) ms\n".format(mean(dtax), std(dtax), min(dtax), max(dtax))
                  text+="p_fail = {} +- {}(min: {}, max: {})\n".format(mean(p_fail), std(p_fail), min(p_fail), max(p_fail))
                  text+='-'*50 + '\n'
                  
                  print(text, file=f)
                    
            if len(GABA_idc)> 0:
                use = SynP[1, GABA_idc]
                tau_rec = SynP[2, GABA_idc]
                tau_fac = SynP[3, GABA_idc]
                g_max = SynP[4, GABA_idc]
                dtax = SynP[5, GABA_idc]
                p_fail = SynP[6, GABA_idc]
                
                if not bins:
                    bins = int(round(np.sqrt(len(GABA_idc)), 0))
                    
                fig, ax = subplots()
                ax.hist(use, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('U')
                fig.savefig('{}/SynapsesParams/GABA_U_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(tau_rec, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('tau_rec (ms)')
                fig.savefig('{}/SynapsesParams/GABA_tau_rec_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(tau_fac, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('tau_fac (ms)')
                fig.savefig('{}/SynapsesParams/GABA_tau_fac_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(g_max, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('g_max (nS)')
                fig.savefig('{}/SynapsesParams/GABA_g_max_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(dtax, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('tau_D (ms)')
                fig.savefig('{}/SynapsesParams/GABA_tau_D_hist_{}.png'.format(sim_dir, curr_idc))
                            
                fig, ax = subplots()
                ax.hist(p_fail, bins=bins)
                ax.set_ylabel('Number of synapses')
                ax.set_xlabel('p_fail')
                fig.savefig('{}/SynapsesParams/GABA_p_fail_hist_{}.png'.format(sim_dir, curr_idc))
                            
                    
                
    
          
                with open(NPfile, 'a') as f:
                    
                    text = "GABA synapses\n\n"
                    
                    text+="Number of synapses: {}\n".format(len(GABA_idc))
                    text+="Mean pCon: {}\n".format(len(GABA_idc)/(Ntargetxsource))
                    
                    text+="Mean +- std (min, max)\n\n"
        
                    text+="U = {} +- {} (min: {}, max: {})\n".format(mean(use), std(use), min(use), max(use))
                    text+="tau_rec = {} +- {} (min: {}, max: {}) ms\n".format(mean(tau_rec), std(tau_rec), min(tau_rec), max(tau_rec))
                    text+="tau_fac = {} +- {} (min: {}, max: {}) ms\n\n".format(mean(tau_fac), std(tau_fac), min(tau_fac), max(tau_fac))
                    
                    text+="g_max = {} +- {} (min: {}, max: {}) nS\n".format(mean(g_max), std(g_max), min(g_max), max(g_max))
                    text+="delay = {} +- {} (min: {}, max: {}) ms\n".format(mean(dtax), std(dtax), min(dtax), max(dtax))
                    text+="p_fail = {} +- {}(min: {}, max: {})\n".format(mean(p_fail), std(p_fail), min(p_fail), max(p_fail))
                    text+='-'*50 + '\n'
    
                    print(text, file=f)
                    
            with open(NPfile, 'a') as f:
                f.write('-'*35)
      
                
    def neuron_t_test(self, pop1, pop2, sim_dir):
        
        if not os.path.isdir('{}/Neuronttest'.format(sim_dir)):
            os.mkdir('{}/Neuronttest'.format(sim_dir))
        
            
         
            
        curr_idc = 0
        NPfile = '{}/Neuronttest/Neuronttest_{}.txt'.format(sim_dir, curr_idc)
        while os.path.isfile(NPfile):
            curr_idc += 1
            NPfile = '{}/Neuronttest/Neuronttest_{}.txt'.format(sim_dir, curr_idc)
            
        
        C1 = self.NeuPar[0, pop1]
        g_L1 = self.NeuPar[1, pop1]
        E_L1 = self.NeuPar[2, pop1]
        delta_T1 = self.NeuPar[3, pop1]
        V_up1 = self.NeuPar[4, pop1]
        tau_w1 = self.NeuPar[5, pop1]
        b1= self.NeuPar[7, pop1]
        V_r1 = self.NeuPar[8, pop1]
        V_T1= self.NeuPar[9, pop1]
        I_ref1 = self.NeuPar[10, pop1]
        V_dep1 = self.NeuPar[11, pop1]
        tau_m1 = self.NeuPar[0, pop1]/self.NeuPar[1, pop1]        
        I_SN1 = self.I_SN[pop1]/pA
        
        C2 = self.NeuPar[0, pop2]
        g_L2 = self.NeuPar[1, pop2]
        E_L2 = self.NeuPar[2, pop2]
        delta_T2 = self.NeuPar[3, pop2]
        V_up2 = self.NeuPar[4, pop2]
        tau_w2 = self.NeuPar[5, pop2]
        b2= self.NeuPar[7, pop2]
        V_r2 = self.NeuPar[8, pop2]
        V_T2= self.NeuPar[9, pop2]
        I_ref2 = self.NeuPar[10, pop2]
        V_dep2 = self.NeuPar[11, pop2]
        tau_m2 = self.NeuPar[0, pop2]/self.NeuPar[1, pop2]        
        I_SN2 = self.I_SN[pop2]/pA
        
        Ctt0, Ctt1 = ttest(C1, C2)
        g_Ltt0, g_Ltt1 = ttest(g_L1, g_L2)
        E_Ltt0, E_Ltt1 = ttest(E_L1, E_L2)
        delta_Ttt0, delta_Ttt1 = ttest(delta_T1, delta_T2)
        V_uptt0, V_uptt1 = ttest(V_up1, V_up2)
        tau_wtt0, tau_wtt1 = ttest(tau_w1, tau_w2)
        btt0, btt1 = ttest(b1, b2)
        V_rtt0, V_rtt1 = ttest(V_r1, V_r2)
        V_Ttt0, V_Ttt1 = ttest(V_T1, V_T2)
        I_reftt0, I_reftt1 = ttest(I_ref1, I_ref2)
        V_deptt0, V_deptt1 = ttest(V_dep1, V_dep2)
        tau_mtt0, tau_mtt1 = ttest(tau_m1, tau_m2)
        I_SNtt0, I_SNtt1 = ttest(I_SN1, I_SN2)
       

        text = ""
        with open(NPfile, 'a') as f:
            text += "t-test of membrane parameters\n\n"
            
            text+= 'pop1: {}\n\n pop2: {}\n\n'.format(pop1, pop2)
            
            text+="C\nstatistic: {}, p-value: {}\n\n".format(round(Ctt0, 3), round(Ctt1, 3))
            text+="g_L\nstatistic: {}, p-value: {}\n\n".format(round(g_Ltt0, 3), round(g_Ltt1, 3))
            text+="E_L\nstatistic: {}, p-value: {}\n\n".format(round(E_Ltt0, 3), round(E_Ltt1, 3))
            text+="delta_T\nstatistic: {}, p-value: {}\n\n".format(round(delta_Ttt0, 3), round(delta_Ttt1, 3))
            text+="V_up\nstatistic: {}, p-value: {}\n\n".format(round(V_uptt0, 3), round(V_uptt1, 3))
            text+="tau_w\nstatistic: {}, p-value: {}\n\n".format(round(tau_wtt0, 3), round(tau_wtt1, 3))
            text+="b\nstatistic: {}, p-value: {}\n\n".format(round(btt0, 3), round(btt1, 3))
            text+="V_r\nstatistic: {}, p-value: {}\n\n".format(round(V_rtt0, 3), round(V_rtt1, 3))
            text+="V_T\nstatistic: {}, p-value: {}\n\n".format(round(V_Ttt0, 3), round(V_Ttt1, 3))
            text+="I_ref\nstatistic: {}, p-value: {}\n\n".format(round(I_reftt0, 3), round(I_reftt1, 3))
            text+="V_dep\nstatistic: {}, p-value: {}\n\n".format(round(V_deptt0, 3), round(V_deptt1, 3))
            text+="tau_m\nstatistic: {}, p-value: {}\n\n".format(round(tau_mtt0, 3), round(tau_mtt1, 3))
            text+="I_SN\nstatistic: {}, p-value: {}\n\n".format(round(I_SNtt0, 3), round(I_SNtt1, 3))
                
            print(text, file=f)

    def neuron_mw_test(self, pop1, pop2, sim_dir):
        
        if not os.path.isdir('{}/Neuronmwtest'.format(sim_dir)):
            os.mkdir('{}/Neuronmwtest'.format(sim_dir))
        
            
         
            
        curr_idc = 0
        NPfile = '{}/Neuronmwtest/Neuronmwtest_{}.txt'.format(sim_dir, curr_idc)
        while os.path.isfile(NPfile):
            curr_idc += 1
            NPfile = '{}/Neuronmwtest/Neuronmwtest_{}.txt'.format(sim_dir, curr_idc)
            
        
        C1 = self.NeuPar[0, pop1]
        g_L1 = self.NeuPar[1, pop1]
        E_L1 = self.NeuPar[2, pop1]
        delta_T1 = self.NeuPar[3, pop1]
        V_up1 = self.NeuPar[4, pop1]
        tau_w1 = self.NeuPar[5, pop1]
        b1= self.NeuPar[7, pop1]
        V_r1 = self.NeuPar[8, pop1]
        V_T1= self.NeuPar[9, pop1]
        I_ref1 = self.NeuPar[10, pop1]
        V_dep1 = self.NeuPar[11, pop1]
        tau_m1 = self.NeuPar[0, pop1]/self.NeuPar[1, pop1]        
        I_SN1 = self.I_SN[pop1]/pA
        
        C2 = self.NeuPar[0, pop2]
        g_L2 = self.NeuPar[1, pop2]
        E_L2 = self.NeuPar[2, pop2]
        delta_T2 = self.NeuPar[3, pop2]
        V_up2 = self.NeuPar[4, pop2]
        tau_w2 = self.NeuPar[5, pop2]
        b2= self.NeuPar[7, pop2]
        V_r2 = self.NeuPar[8, pop2]
        V_T2= self.NeuPar[9, pop2]
        I_ref2 = self.NeuPar[10, pop2]
        V_dep2 = self.NeuPar[11, pop2]
        tau_m2 = self.NeuPar[0, pop2]/self.NeuPar[1, pop2]        
        I_SN2 = self.I_SN[pop2]/pA
        
        Cmw0, Cmw1 = mwtest(C1, C2, alternative='two-sided')
        g_Lmw0, g_Lmw1 = mwtest(g_L1, g_L2, alternative='two-sided')
        E_Lmw0, E_Lmw1 = mwtest(E_L1, E_L2, alternative='two-sided')
        delta_mwt0, delta_mwt1 = mwtest(delta_T1, delta_T2, alternative='two-sided')
        V_upmw0, V_upmw1 = mwtest(V_up1, V_up2, alternative='two-sided')
        tau_wmw0, tau_wmw1 = mwtest(tau_w1, tau_w2, alternative='two-sided')
        bmw0, bmw1 = mwtest(b1, b2, alternative='two-sided')
        V_rmw0, V_rmw1 = mwtest(V_r1, V_r2, alternative='two-sided')
        V_mwt0, V_mwt1 = mwtest(V_T1, V_T2, alternative='two-sided')
        I_refmw0, I_refmw1 = mwtest(I_ref1, I_ref2, alternative='two-sided')
        V_depmw0, V_depmw1 = mwtest(V_dep1, V_dep2, alternative='two-sided')
        tau_mmw0, tau_mmw1 = mwtest(tau_m1, tau_m2, alternative='two-sided')
        I_SNmw0, I_SNmw1 = mwtest(I_SN1, I_SN2, alternative='two-sided')
       

        text = ""
        with open(NPfile, 'a') as f:
            text += "Mann-Whitney U test of membrane parameters\n\n"
            
            text+= 'pop1: {}\n\n pop2: {}\n\n'.format(pop1, pop2)
            
            text+="C\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(Cmw0, 3), round(Cmw1, 3), round(np.mean(C1), 3), round(np.std(C1), 3), round(np.mean(C2), 3), round(np.std(C2), 3))
            text+="g_L\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(g_Lmw0, 3), round(g_Lmw1, 3), round(np.mean(g_L1), 3), round(np.std(g_L1), 3), round(np.mean(g_L2), 3), round(np.std(g_L2), 3))
            text+="E_L\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(E_Lmw0, 3), round(E_Lmw1, 3), round(np.mean(E_L1), 3), round(np.std(E_L1), 3), round(np.mean(E_L2), 3), round(np.std(E_L2), 3))
            text+="delta_T\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(delta_mwt0, 3), round(delta_mwt1, 3), round(np.mean(delta_T1), 3), round(np.std(delta_T1), 3), round(np.mean(delta_T2), 3), round(np.std(delta_T2), 3))
            text+="V_up\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(V_upmw0, 3), round(V_upmw1, 3), round(np.mean(V_up1), 3), round(np.std(V_up1), 3), round(np.mean(V_up2), 3), round(np.std(V_up2), 3))
            text+="tau_w\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(tau_wmw0, 3), round(tau_wmw1, 3), round(np.mean(tau_w1), 3), round(np.std(tau_w1), 3), round(np.mean(tau_w2), 3), round(np.std(tau_w2), 3))
            text+="b\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(bmw0, 3), round(bmw1, 3), round(np.mean(b1), 3), round(np.std(b1), 3), round(np.mean(b2), 3), round(np.std(b2), 3))
            text+="V_r\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(V_rmw0, 3), round(V_rmw1, 3), round(np.mean(V_r1), 3), round(np.std(V_r1), 3), round(np.mean(V_r2), 3), round(np.std(V_r2), 3))
            text+="V_T\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(V_mwt0, 3), round(V_mwt1, 3), round(np.mean(V_T1), 3), round(np.std(V_T1), 3), round(np.mean(V_T2), 3), round(np.std(V_T2), 3))
            text+="I_ref\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(I_refmw0, 3), round(I_refmw1, 3), round(np.mean(I_ref1), 3), round(np.std(I_ref1), 3), round(np.mean(I_ref2), 3), round(np.std(I_ref2), 3))
            text+="V_dep\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(V_depmw0, 3), round(V_depmw1, 3), round(np.mean(V_dep1), 3), round(np.std(V_dep1), 3), round(np.mean(V_dep2), 3), round(np.std(V_dep2), 3))
            text+="tau_m\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(tau_mmw0, 3), round(tau_mmw1, 3), round(np.mean(tau_m1), 3), round(np.std(tau_m1), 3), round(np.mean(tau_m2), 3), round(np.std(tau_m2), 3))
            text+="I_SN\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(I_SNmw0, 3), round(I_SNmw1, 3), round(np.mean(I_SN1), 3), round(np.std(I_SN1), 3), round(np.mean(I_SN2), 3), round(np.std(I_SN2), 3))
                
            print(text, file=f)
        
    def synapses_t_test(self, pop1, pop2, sim_dir, target=True):
      
        if not os.path.isdir('{}/Synapsesttest'.format(sim_dir)):
            os.mkdir('{}/Synapsesttest'.format(sim_dir))
        
            
       
        curr_idc = 0
        NPfile = '{}/Synapsesttest/Synapsesttest_{}.txt'.format(sim_dir, curr_idc)
        
        while os.path.isfile(NPfile):
            curr_idc += 1
            NPfile = '{}/Synapsesttest/Synapsesttest_{}.txt'.format(sim_dir, curr_idc)
            
        if target:
            syn_idc1 = np.asarray(self.group_name_to_synapses_idc([[pop1, 'all']], neuron_index=True))          
            syn_idc2 = np.asarray(self.group_name_to_synapses_idc([[pop2, 'all']], neuron_index=True))  
        else:
            syn_idc1 = np.asarray(self.group_name_to_synapses_idc([['all', pop1]], neuron_index=True))          
            syn_idc2 = np.asarray(self.group_name_to_synapses_idc([['all', pop2]], neuron_index=True))

            
        Ntarget1 = len(pop1)
        Nsource1 = self.NeuPar.shape[1]
        Ntargetxsource1 = Ntarget1*Nsource1
        
        SynP1 = self.SynPar[:, syn_idc1]
        stypes1 = SynP1[0, :]
        
        AMPA_idc1 = np.where(stypes1==1)[0]
        NMDA_idc1 = np.where(stypes1==3)[0]
        GABA_idc1 = np.where(stypes1==2)[0]

        
        Ntarget2 = len(pop2)
        Nsource2 = self.NeuPar.shape[1]
        Ntargetxsource2 = Ntarget2*Nsource2
        
        SynP2 = self.SynPar[:, syn_idc2]
        stypes2 = SynP2[0, :]
        
        AMPA_idc2 = np.where(stypes2==1)[0]
        NMDA_idc2 = np.where(stypes2==3)[0]
        GABA_idc2 = np.where(stypes2==2)[0]
        
        with open(NPfile, 'a') as f:
            text = "synapses t-test\n\n"
            text += "pop1: {}\n\npop2: {}\n\n".format(pop1, pop2)
            print(text, file=f)
            
        if len(AMPA_idc1)*len(AMPA_idc2) > 0:
            
            use1 = SynP1[1, AMPA_idc1]
            tau_rec1 = SynP1[2, AMPA_idc1]
            tau_fac1 = SynP1[3, AMPA_idc1]
            g_max1 = SynP1[4, AMPA_idc1]
            dtax1 = SynP1[5, AMPA_idc1]
            p_fail1 = SynP1[6, AMPA_idc1]
            
            use2 = SynP2[1, AMPA_idc2]
            tau_rec2 = SynP2[2, AMPA_idc2]
            tau_fac2 = SynP2[3, AMPA_idc2]
            g_max2 = SynP2[4, AMPA_idc2]
            dtax2 = SynP2[5, AMPA_idc2]
            p_fail2 = SynP2[6, AMPA_idc2]
            
            usett1, usett2 = ttest(use1, use2)
            tau_rectt1, tau_rectt2 = ttest(tau_rec1, tau_rec2)
            tau_factt1, tau_factt2 = ttest(tau_fac1, tau_fac2)
            g_maxtt1, g_maxtt2 = ttest(g_max1, g_max2)
            dtaxtt1, dtaxtt2 = ttest(dtax1, dtax2)
            p_failtt1, p_failtt2 = ttest(p_fail1, p_fail2)
            
            
            with open(NPfile, 'a') as f:
                
              text = "AMPA synapses\n\n"
              
              text += "U\nstatistic: {}, p-value: {}\n\n".format(round(usett1, 3), round(usett2, 3))
              text += "tau_rec\nstatistic: {}, p-value: {}\n\n".format(round(tau_rectt1, 3), round(tau_rectt2, 3))
              text += "tau_fac\nstatistic: {}, p-value: {}\n\n".format(round(tau_factt1, 3), round(tau_factt2, 3))
              text += "g_max\nstatistic: {}, p-value: {}\n\n".format(round(g_maxtt1, 3), round(g_maxtt2, 3))
              text += "dtax\nstatistic: {}, p-value: {}\n\n".format(round(dtaxtt1, 3), round(dtaxtt2, 3))
              text += "p_fail\nstatistic: {}, p-value: {}\n\n".format(round(p_failtt1, 3), round(p_failtt2, 3))
              
              print(text, file=f)
              
        if len(NMDA_idc1)*len(NMDA_idc2) > 0:
            use1 = SynP1[1, NMDA_idc1]
            tau_rec1 = SynP1[2, NMDA_idc1]
            tau_fac1 = SynP1[3, NMDA_idc1]
            g_max1 = SynP1[4, NMDA_idc1]
            dtax1 = SynP1[5, NMDA_idc1]
            p_fail1 = SynP1[6, NMDA_idc1]
            
            use2 = SynP2[1, NMDA_idc2]
            tau_rec2 = SynP2[2, NMDA_idc2]
            tau_fac2 = SynP2[3, NMDA_idc2]
            g_max2 = SynP2[4, NMDA_idc2]
            dtax2 = SynP2[5, NMDA_idc2]
            p_fail2 = SynP2[6, NMDA_idc2]
            
            usett1, usett2 = ttest(use1, use2)
            tau_rectt1, tau_rectt2 = ttest(tau_rec1, tau_rec2)
            tau_factt1, tau_factt2 = ttest(tau_fac1, tau_fac2)
            g_maxtt1, g_maxtt2 = ttest(g_max1, g_max2)
            dtaxtt1, dtaxtt2 = ttest(dtax1, dtax2)
            p_failtt1, p_failtt2 = ttest(p_fail1, p_fail2)
            
            
            with open(NPfile, 'a') as f:
                
              text = "NMDA synapses\n\n"
              
              text += "U\nstatistic: {}, p-value: {}\n\n".format(round(usett1, 3), round(usett2, 3))
              text += "tau_rec\nstatistic: {}, p-value: {}\n\n".format(round(tau_rectt1, 3), round(tau_rectt2, 3))
              text += "tau_fac\nstatistic: {}, p-value: {}\n\n".format(round(tau_factt1, 3), round(tau_factt2, 3))
              text += "g_max\nstatistic: {}, p-value: {}\n\n".format(round(g_maxtt1, 3), round(g_maxtt2, 3))
              text += "dtax\nstatistic: {}, p-value: {}\n\n".format(round(dtaxtt1, 3), round(dtaxtt2, 3))
              text += "p_fail\nstatistic: {}, p-value: {}\n\n".format(round(p_failtt1, 3), round(p_failtt2, 3))
              
              print(text, file=f)
            
            
            
        if len(GABA_idc1)*len(GABA_idc2) > 0:
            use1 = SynP1[1, GABA_idc1]
            tau_rec1 = SynP1[2, GABA_idc1]
            tau_fac1 = SynP1[3, GABA_idc1]
            g_max1 = SynP1[4, GABA_idc1]
            dtax1 = SynP1[5, GABA_idc1]
            p_fail1 = SynP1[6, GABA_idc1]
            
            use2 = SynP2[1, GABA_idc2]
            tau_rec2 = SynP2[2, GABA_idc2]
            tau_fac2 = SynP2[3, GABA_idc2]
            g_max2 = SynP2[4, GABA_idc2]
            dtax2 = SynP2[5, GABA_idc2]
            p_fail2 = SynP2[6, GABA_idc2]
            
            usett1, usett2 = ttest(use1, use2)
            tau_rectt1, tau_rectt2 = ttest(tau_rec1, tau_rec2)
            tau_factt1, tau_factt2 = ttest(tau_fac1, tau_fac2)
            g_maxtt1, g_maxtt2 = ttest(g_max1, g_max2)
            dtaxtt1, dtaxtt2 = ttest(dtax1, dtax2)
            p_failtt1, p_failtt2 = ttest(p_fail1, p_fail2)
            
            
            with open(NPfile, 'a') as f:
                
              text = "GABA synapses\n\n"
              
              text += "U\nstatistic: {}, p-value: {}\n\n".format(round(usett1, 3), round(usett2, 3))
              text += "tau_rec\nstatistic: {}, p-value: {}\n\n".format(round(tau_rectt1, 3), round(tau_rectt2, 3))
              text += "tau_fac\nstatistic: {}, p-value: {}\n\n".format(round(tau_factt1, 3), round(tau_factt2, 3))
              text += "g_max\nstatistic: {}, p-value: {}\n\n".format(round(g_maxtt1, 3), round(g_maxtt2, 3))
              text += "dtax\nstatistic: {}, p-value: {}\n\n".format(round(dtaxtt1, 3), round(dtaxtt2, 3))
              text += "p_fail\nstatistic: {}, p-value: {}\n\n".format(round(p_failtt1, 3), round(p_failtt2, 3))
              
              print(text, file=f)
        
    def synapses_mw_test(self, pop1, pop2, sim_dir, target=True):
      
        if not os.path.isdir('{}/Synapsesmwtest'.format(sim_dir)):
            os.mkdir('{}/Synapsesmwtest'.format(sim_dir))
        
            
       
        curr_idc = 0
        NPfile = '{}/Synapsesmwtest/Synapsesmwtest_{}.txt'.format(sim_dir, curr_idc)
        
        while os.path.isfile(NPfile):
            curr_idc += 1
            NPfile = '{}/Synapsesmwtest/Synapsesmwtest_{}.txt'.format(sim_dir, curr_idc)
            
        if target:
            syn_idc1 = np.asarray(self.group_name_to_synapses_idc([[pop1, 'all']], neuron_index=True))          
            syn_idc2 = np.asarray(self.group_name_to_synapses_idc([[pop2, 'all']], neuron_index=True))  
        else:
            syn_idc1 = np.asarray(self.group_name_to_synapses_idc([['all', pop1]], neuron_index=True))          
            syn_idc2 = np.asarray(self.group_name_to_synapses_idc([['all', pop2]], neuron_index=True))

            
        Ntarget1 = len(pop1)
        Nsource1 = self.NeuPar.shape[1]
        Ntargetxsource1 = Ntarget1*Nsource1
        
        SynP1 = self.SynPar[:, syn_idc1]
        stypes1 = SynP1[0, :]
        
        AMPA_idc1 = np.where(stypes1==1)[0]
        NMDA_idc1 = np.where(stypes1==3)[0]
        GABA_idc1 = np.where(stypes1==2)[0]

        
        Ntarget2 = len(pop2)
        Nsource2 = self.NeuPar.shape[1]
        Ntargetxsource2 = Ntarget2*Nsource2
        
        SynP2 = self.SynPar[:, syn_idc2]
        stypes2 = SynP2[0, :]
        
        AMPA_idc2 = np.where(stypes2==1)[0]
        NMDA_idc2 = np.where(stypes2==3)[0]
        GABA_idc2 = np.where(stypes2==2)[0]
        
        with open(NPfile, 'a') as f:
            text = "synapses Mann-Whitney U-test\n\n"
            text += "pop1: {}\n\npop2: {}\n\n".format(pop1, pop2)
            print(text, file=f)
            
        if len(AMPA_idc1)*len(AMPA_idc2) > 0:
            
            use1 = SynP1[1, AMPA_idc1]
            tau_rec1 = SynP1[2, AMPA_idc1]
            tau_fac1 = SynP1[3, AMPA_idc1]
            g_max1 = SynP1[4, AMPA_idc1]
            dtax1 = SynP1[5, AMPA_idc1]

            use2 = SynP2[1, AMPA_idc2]
            tau_rec2 = SynP2[2, AMPA_idc2]
            tau_fac2 = SynP2[3, AMPA_idc2]
            g_max2 = SynP2[4, AMPA_idc2]
            dtax2 = SynP2[5, AMPA_idc2]

            
            usemw1, usemw2 = mwtest(use1, use2, alternative='two-sided')
            tau_recmw1, tau_recmw2 = mwtest(tau_rec1, tau_rec2, alternative='two-sided')
            tau_facmw1, tau_facmw2 = mwtest(tau_fac1, tau_fac2, alternative='two-sided')
            g_maxmw1, g_maxmw2 = mwtest(g_max1, g_max2, alternative='two-sided')
            dtaxmw1, dtaxmw2 = mwtest(dtax1, dtax2, alternative='two-sided')
     
            
            
            with open(NPfile, 'a') as f:
                
                text = "AMPA synapses\n\n"
                
                text += "U\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(usemw1, 3), round(usemw2, 3), round(np.mean(use1), 3), round(np.std(use1), 3), round(np.mean(use2), 3), round(np.std(use2), 3))
                text += "tau_rec\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(tau_recmw1, 3), round(tau_recmw2, 3), round(np.mean(tau_rec1), 3), round(np.std(tau_rec1), 3), round(np.mean(tau_rec2), 3), round(np.std(tau_rec2), 3))
                text += "tau_fac\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(tau_facmw1, 3), round(tau_facmw2, 3), round(np.mean(tau_fac1), 3), round(np.std(tau_fac1), 3), round(np.mean(tau_fac2), 3), round(np.std(tau_fac2), 3))
                text += "g_max\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(g_maxmw1, 3), round(g_maxmw2, 3), round(np.mean(g_max1), 3), round(np.std(g_max1), 3), round(np.mean(g_max2), 3), round(np.std(g_max2), 3))
                text += "dtax\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(dtaxmw1, 3), round(dtaxmw2, 3), round(np.mean(dtax1), 3), round(np.std(dtax1), 3), round(np.mean(dtax2), 3), round(np.std(dtax2), 3))
         
                
                print(text, file=f)
              
        if len(NMDA_idc1)*len(NMDA_idc2) > 0:
            use1 = SynP1[1, NMDA_idc1]
            tau_rec1 = SynP1[2, NMDA_idc1]
            tau_fac1 = SynP1[3, NMDA_idc1]
            g_max1 = SynP1[4, NMDA_idc1]
            dtax1 = SynP1[5, NMDA_idc1]

            
            use2 = SynP2[1, NMDA_idc2]
            tau_rec2 = SynP2[2, NMDA_idc2]
            tau_fac2 = SynP2[3, NMDA_idc2]
            g_max2 = SynP2[4, NMDA_idc2]
            dtax2 = SynP2[5, NMDA_idc2]

            
            usemw1, usemw2 = mwtest(use1, use2, alternative='two-sided')
            tau_recmw1, tau_recmw2 = mwtest(tau_rec1, tau_rec2, alternative='two-sided')
            tau_facmw1, tau_facmw2 = mwtest(tau_fac1, tau_fac2, alternative='two-sided')
            g_maxmw1, g_maxmw2 = mwtest(g_max1, g_max2, alternative='two-sided')
            dtaxmw1, dtaxmw2 = mwtest(dtax1, dtax2, alternative='two-sided')

            
            
            with open(NPfile, 'a') as f:
                
                text = "NMDA synapses\n\n"
                
                text += "U\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(usemw1, 3), round(usemw2, 3), round(np.mean(use1), 3), round(np.std(use1), 3), round(np.mean(use2), 3), round(np.std(use2), 3))
                text += "tau_rec\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(tau_recmw1, 3), round(tau_recmw2, 3), round(np.mean(tau_rec1), 3), round(np.std(tau_rec1), 3), round(np.mean(tau_rec2), 3), round(np.std(tau_rec2), 3))
                text += "tau_fac\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(tau_facmw1, 3), round(tau_facmw2, 3), round(np.mean(tau_fac1), 3), round(np.std(tau_fac1), 3), round(np.mean(tau_fac2), 3), round(np.std(tau_fac2), 3))
                text += "g_max\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(g_maxmw1, 3), round(g_maxmw2, 3), round(np.mean(g_max1), 3), round(np.std(g_max1), 3), round(np.mean(g_max2), 3), round(np.std(g_max2), 3))
                text += "dtax\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(dtaxmw1, 3), round(dtaxmw2, 3), round(np.mean(dtax1), 3), round(np.std(dtax1), 3), round(np.mean(dtax2), 3), round(np.std(dtax2), 3))
         
              
                print(text, file=f)
            
            
            
        if len(GABA_idc1)*len(GABA_idc2) > 0:
            use1 = SynP1[1, GABA_idc1]
            tau_rec1 = SynP1[2, GABA_idc1]
            tau_fac1 = SynP1[3, GABA_idc1]
            g_max1 = SynP1[4, GABA_idc1]
            dtax1 = SynP1[5, GABA_idc1]

            
            use2 = SynP2[1, GABA_idc2]
            tau_rec2 = SynP2[2, GABA_idc2]
            tau_fac2 = SynP2[3, GABA_idc2]
            g_max2 = SynP2[4, GABA_idc2]
            dtax2 = SynP2[5, GABA_idc2]

            
            usemw1, usemw2 = mwtest(use1, use2, alternative='two-sided')
            tau_recmw1, tau_recmw2 = mwtest(tau_rec1, tau_rec2, alternative='two-sided')
            tau_facmw1, tau_facmw2 = mwtest(tau_fac1, tau_fac2, alternative='two-sided')
            g_maxmw1, g_maxmw2 = mwtest(g_max1, g_max2, alternative='two-sided')
            dtaxmw1, dtaxmw2 = mwtest(dtax1, dtax2, alternative='two-sided')

            
            
            with open(NPfile, 'a') as f:
                 
                text = "GABA synapses\n\n"
              
                text += "U\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(usemw1, 3), round(usemw2, 3), round(np.mean(use1), 3), round(np.std(use1), 3), round(np.mean(use2), 3), round(np.std(use2), 3))
                text += "tau_rec\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(tau_recmw1, 3), round(tau_recmw2, 3), round(np.mean(tau_rec1), 3), round(np.std(tau_rec1), 3), round(np.mean(tau_rec2), 3), round(np.std(tau_rec2), 3))
                text += "tau_fac\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(tau_facmw1, 3), round(tau_facmw2, 3), round(np.mean(tau_fac1), 3), round(np.std(tau_fac1), 3), round(np.mean(tau_fac2), 3), round(np.std(tau_fac2), 3))
                text += "g_max\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(g_maxmw1, 3), round(g_maxmw2, 3), round(np.mean(g_max1), 3), round(np.std(g_max1), 3), round(np.mean(g_max2), 3), round(np.std(g_max2), 3))
                text += "dtax\nstatistic: {}, p-value: {}\npop1: {} +- {} (mean +- std)\npop2: {} +- {} (mean +- std)\n\n".format(round(dtaxmw1, 3), round(dtaxmw2, 3), round(np.mean(dtax1), 3), round(np.std(dtax1), 3), round(np.mean(dtax2), 3), round(np.std(dtax2), 3))
         
              
                print(text, file=f)
        
        
    def connection_chi(self, pop1, pop2, sim_dir, synapses_type='all'):
        
        if not os.path.isdir('{}/pCon_chitest'.format(sim_dir)):
            os.mkdir('{}/pCon_chitest'.format(sim_dir))
        
            
       
        curr_idc = 0
        NPfile = '{}/pCon_chitest/pCon_chitest_{}.txt'.format(sim_dir, curr_idc)
        
        while os.path.isfile(NPfile):
            curr_idc += 1
            NPfile = '{}/pCon_chitest/pCon_chitest_{}.txt'.format(sim_dir, curr_idc)
            
        
        if synapses_type=='all':
            pCon1 = self.find_pCon([[pop1, 'all']], neuron_index=True)
            pCon2 = self.find_pCon([[pop2, 'all']], neuron_index=True)
            
            Ntotal = self.NeuPar.shape[1]
            
            Npop1 = int(pCon1 * Ntotal)
            Npop1comp = int((1-pCon1) * Ntotal)
            Npop2 = int(pCon2 * Ntotal)
            Npop2comp = int((1-pCon2)*Ntotal)
            
        elif synapses_type=='excitatory':
            pCon1 = self.find_exc_pCon([[pop1, 'all']], neuron_index=True)
            pCon2 = self.find_exc_pCon([[pop2, 'all']], neuron_index=True)
            
            excidc=[]
            for col in range(len(self.group_distr[0])):
                excidc.extend(list(self.group_name_to_neuron_idc([['PC', col]])))
            
            Ntotal = len(excidc)
            
            Npop1 = int(pCon1 * Ntotal)
            Npop1comp = int((1-pCon1) * Ntotal)
            Npop2 = int(pCon2 * Ntotal)
            Npop2comp = int((1-pCon2)*Ntotal)
            
        elif synapses_type=='inhibitory':
            pCon1 = self.find_inh_pCon([[pop1, 'all']], neuron_index=True)
            pCon2 = self.find_inh_pCon([[pop2, 'all']], neuron_index=True)
            
            inhidc=[]
            for col in range(len(self.group_distr[0])):
                inhidc.extend(list(self.group_name_to_neuron_idc([['IN', col]])))
            
            Ntotal = len(inhidc)
            
            Npop1 = int(pCon1 * Ntotal)
            Npop1comp = int((1-pCon1) * Ntotal)
            Npop2 = int(pCon2 * Ntotal)
            Npop2comp = int((1-pCon2)*Ntotal)
            
        chi_, p_value, _1, _2 = chi2([[Npop1, Npop1comp], [Npop2, Npop2comp]])
            
        with open(NPfile, 'a') as f:
            print('Squared chi test for connectivity - {} synapses\n'.format(synapses_type), file=f)
            
            print('pCon pop1', pCon1, file=f)
            print('pCon pop2', pCon2, end='\n\n', file=f)
            
            print('Chi:', chi_, file=f)
            print('p-value', p_value, file=f)
            
        
        return chi_, p_value, pCon1, pCon2
            
            
            
            
            
            # excidc = []
            # for col in range(len(self.group_distr[0])):
            #     excidc.extend(self.group_name_to_neuron_idc([['PC', 0]]))
    
        
    
    def ISI_analysis(self, analysis, sim_dir):

        print('REPORT: Starting ISI analysis\n')
        if not os.path.isdir('{}/ISI_analysis'.format(sim_dir)):
            os.mkdir('{}/ISI_analysis'.format(sim_dir))      
       
        zerolagCC_list = []
        ISImean_list = []
        ISIstd_list = []
        ISICV_list = []
        spiketime_list = []
        return_expecvalue_list = []
        return_var_list = []
        autocorrvalue_list = []
        autocorrt_list = []
        binned_spiketimes_list = []
        
        t_arr = self.spikemonitor.t/ms
     
        for analysand in analysis:  
            
            curr_idc = 0
            
            ISI_file = '{}/ISI_analysis/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
            
            while os.path.isfile(ISI_file):
                curr_idc+= 1
                ISI_file = '{}/ISI_analysis/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
             
            print('Analysis {}'.format(curr_idc))
            
            group_n_stripes = analysand['Group']
            start_time = analysand['Start']
            stop_time = analysand['Stop']
            time_bin = analysand['Time bin']
            min_sp_num = analysand['Minimum spike number']
                
            neurons = self.group_name_to_neuron_idc(group_n_stripes)
            idc_duration = np.where((self.spikemonitor.t/ms>= start_time)&(self.spikemonitor.t/ms<stop_time))[0]
            neurons_duration, counts = np.unique(self.spikemonitor.i[idc_duration], return_counts=True)
            
            spiking = neurons_duration[np.where(counts>= min_sp_num)]
            neurons_idc = np.intersect1d(neurons, spiking)
             
            if len(neurons_idc):
                time_interval = stop_time - start_time
                time_interval_in_s = time_interval/1000
                time_bin_in_s = time_bin/1000
                N_bins = int(round(time_interval/time_bin, 0))
                
                expecvalue_list = []
                var_list = []
                r_list = []
                
                ISI_mean = []
                ISI_std = []
                ISI_CV = []
                sp_arr = []
                
                binned_spiketimes = np.zeros((len(neurons_idc), N_bins))
                
                k = 0
                dez = 0
                totalneurons = len(neurons_idc)
                for i in range(len(neurons_idc)):
                    k+= 1
                    perc = k/totalneurons*100
                    if perc//10 > dez:
                        dez = perc//10
                        print('{}% of CV concluded'.format(int(round(perc, 0))))
                    T_spike = self.spikemonitor.t[self.spikemonitor.i==neurons_idc[i]]/ms
                    T1 = np.floor((T_spike[np.where((T_spike>=start_time)&(T_spike<= stop_time))]-start_time)/time_bin).astype(int)
    
                    count = len(T1)     
                    
                    ISI = np.diff(T_spike)
                    
                    expecvalue = count/N_bins
                    
                    expecvalue_list.append(expecvalue)
                    var_list.append(expecvalue*(1-expecvalue))
    
                    
                    binned_spiketimes[i, T1] = 1
                    sp_arr.append(T_spike)
                    ISI_mean.append(np.mean(ISI))
                    ISI_std.append(np.std(ISI))
                    ISI_CV.append(np.std(ISI)/np.mean(ISI))
                    
                print('CV concluded\n')
                cov_arr = []
                
                k = 0
                dez = 0
                perc = 0
                totalneurons = len(neurons_idc) * (len(neurons_idc) - 1)//2
                for i in range(len(neurons_idc)):
                    for j in range(i):
                        k+= 1
                        perc = k/totalneurons*100
                        if perc//10 > dez:
                            dez = perc//10
                            print('{}% of cross-correlation concluded'.format(int(round(perc, 0))))
                        p_xy = binned_spiketimes[i, :]@binned_spiketimes[j, :]/N_bins
                        cov_arr.append((    p_xy - expecvalue_list[i]*expecvalue_list[j])/np.sqrt(var_list[i] * var_list[j]))
                        
                print('Cross-correlation concluded\n')
                auto_bins = N_bins//2
                auto_cov_arr = np.zeros((len(neurons_idc), auto_bins))
                
                k = 0
                dez = 0
                totalbins = N_bins//2
                
                for t in range(N_bins//2):
                    k+= 1
                    perc = k/totalbins*100
                    if perc//10 > dez:
                        dez = perc//10
                        print('{}% of autocorrelation concluded'.format(int(round(perc, 0))))
                    if t == 0:
                        b1 = binned_spiketimes[:, :]
                        b2 = binned_spiketimes[:, :]
                        _ = b1*b2
                        p_x = np.sum(binned_spiketimes[:, :], axis=1)/N_bins
                        p_y = p_x
                        
                    else:
                        b1 = binned_spiketimes[:, t:]
                        b2 = binned_spiketimes[:, :-t]
                        _ = b1*b2
                        p_x = np.sum(binned_spiketimes[:, :-t], axis=1)/(N_bins-t)
                        p_y = np.sum(binned_spiketimes[:, t:], axis=1)/(N_bins-t)
               
                    p_xy = np.sum(_, axis=1)/(N_bins-2*t)
                    auto_cov_arr[:, t] = (p_xy - p_x * p_y)/np.sqrt(p_x*(1-p_x) * p_y * (1-p_y))
                    
                print('Autocorrelation concluded\n')
                auto_t = np.arange(auto_bins)*time_bin
                
                with open(ISI_file, 'a') as f:
                    print('Analysis:', analysand,end='\n\n', file=f)
                    
                    print('Spiking neurons (n spikes >= {}):'.format(min_sp_num), len(neurons_idc), end='\n\n', file=f)
                    
                    print('ISI mean mean:', np.mean(ISI_mean), file=f)
                    print('ISI mean std:', np.std(ISI_mean), end='\n\n', file=f)
                    
                    print('ISI std mean:', np.mean(ISI_std), file=f)
                    print('ISI std std:', np.std(ISI_std), end='\n\n', file=f)
                    
                    print('ISI CV mean:', np.mean(ISI_CV), file=f)
                    print('ISI CV std:', np.std(ISI_CV), end='\n\n', file=f)
                    
                    print('CC(0) mean:', np.mean(cov_arr), file=f)
                    print('CC(0) std:', np.std(cov_arr), end='\n\n', file=f)
                 
                zerolagCC_list.append(cov_arr)
                ISImean_list.append(ISI_mean)
                ISIstd_list.append(ISI_std)
                ISICV_list.append(ISI_CV)
                spiketime_list.append(sp_arr)
                return_expecvalue_list.append(expecvalue_list)
                return_var_list .append(var_list)
                autocorrvalue_list.append(auto_cov_arr)
                autocorrt_list.append(auto_t)
                binned_spiketimes_list.append(binned_spiketimes)
                
                
                meanautoC = np.mean(auto_cov_arr, axis=0)
                
    
                fig, ax = subplots()
                ax.plot(auto_t[:len(auto_t)//2], meanautoC[:len(auto_t)//2])
                ax.set_xlabel('Lag (ms)')
                ax.set_ylabel('Autocorrelation')
                fig.savefig('{}/ISI_analysis/meanautoCtrace_{}.png'.format(sim_dir, curr_idc))
                
                
                
                
                auto_t = np.asarray(auto_t) 
                t100ms = auto_t[auto_t<100]
                t100ms = np.where((auto_t<300))[0]
                fig, ax = subplots()
                ax.plot(auto_t[t100ms], meanautoC[t100ms])
                ax.set_xlabel('Lag (ms)')
                ax.set_ylabel('Autocorrelation')
                fig.savefig('{}/ISI_analysis/meanautoCtrace300ms_{}.png'.format(sim_dir, curr_idc))
    
                t10_100ms = np.where((auto_t>=2) & (auto_t<300))[0]
                fig, ax = subplots()
                ax.plot(auto_t[t10_100ms], meanautoC[t10_100ms])
                ax.set_xlabel('Lag (ms)')
                ax.set_ylabel('Autocorrelation')
                fig.savefig('{}/ISI_analysis/meanautoCtrace2_300ms_{}.png'.format(sim_dir, curr_idc))
                
                t10_500ms = np.where((auto_t>=2) & (auto_t<1000))[0]
                fig, ax = subplots()
                ax.plot(auto_t[t10_500ms], meanautoC[t10_500ms])
                ax.set_xlabel('Lag (ms)')
                ax.set_ylabel('Autocorrelation')
                fig.savefig('{}/ISI_analysis/meanautoCtrace2_1000ms_{}.png'.format(sim_dir, curr_idc))
            
            
                bins_lims, dens_lims, dens_arr = build_hist(cov_arr)
                
                fig, ax = subplots()
                ax.hist(cov_arr, bins=int(round(np.sqrt(len(cov_arr)), 0)))
                ax.set_xlabel('Zero-lag cross-correlation')
                ax.set_ylabel('Number of pairs')
                fig.savefig('{}/ISI_analysis/autocorhist_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(cov_arr, bins=bins_lims, density=True)
                ax.set_xlabel('Zero-lag cross-correlation')
                ax.set_ylabel('Density')
                fig.savefig('{}/ISI_analysis/autocorhistdens_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Zero-lag cross-correlation')
                ax.set_ylabel('Density')
                fig.savefig('{}/ISI_analysis/autocorhistplot_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Zero-lag cross-correlation')
                ax.set_ylabel('Density')
                ax.set_xlim(np.mean(cov_arr) - 1*np.std(cov_arr), np.mean(cov_arr) + 1*np.std(cov_arr))
                fig.savefig('{}/ISI_analysis/autocorhistplot1_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Zero-lag cross-correlation')
                ax.set_ylabel('Density')
                ax.set_xlim(np.mean(cov_arr) - 2*np.std(cov_arr), np.mean(cov_arr) + 2*np.std(cov_arr))
                fig.savefig('{}/ISI_analysis/autocorhistplot2_{}.png'.format(sim_dir, curr_idc))
                
                
                bins_lims, dens_lims, dens_arr = build_hist(ISI_CV)

                fig, ax = subplots()
                ax.hist(ISI_CV, bins=int(round(np.sqrt(len(ISI_CV)), 0)))
                ax.set_xlabel('Coefficient of variation')
                ax.set_ylabel('Number of neurons')
                fig.savefig('{}/ISI_analysis/CVhist_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(ISI_CV, bins=bins_lims, density=True)
                ax.set_xlabel('Coefficient of variation')
                ax.set_ylabel('Density')
                fig.savefig('{}/ISI_analysis/CVhistdens_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Coefficient of variation')
                ax.set_ylabel('Density')
                fig.savefig('{}/ISI_analysis/CVhistplot_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Coefficient of variation')
                ax.set_ylabel('Density')
                ax.set_xlim(0, np.mean(ISI_CV))
                fig.savefig('{}/ISI_analysis/CVhistplot0_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Coefficient of variation')
                ax.set_ylabel('Density')
                ax.set_xlim(0, np.mean(ISI_CV) + 1* np.std(ISI_CV))
                fig.savefig('{}/ISI_analysis/CVhistplot1_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Coefficient of variation')
                ax.set_ylabel('Density')
                ax.set_xlim(0, np.mean(ISI_CV) + 2* np.std(ISI_CV))
                fig.savefig('{}/ISI_analysis/CVhistplot2_{}.png'.format(sim_dir, curr_idc))
                
                
                bins_lims, dens_lims, dens_arr = build_hist(ISI_mean)
                
                fig, ax = subplots()
                ax.hist(ISI_mean, bins=int(round(np.sqrt(len(ISI_mean)), 0)))
                ax.set_xlabel('mean ISI (ms)')
                ax.set_ylabel('Number of neurons')
                fig.savefig('{}/ISI_analysis/ISImeanhist_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(ISI_mean, bins=bins_lims, density=True)
                ax.set_xlabel('mean ISI (ms)')
                ax.set_ylabel('Density')
                fig.savefig('{}/ISI_analysis/ISImeanhistdens_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('mean ISI (ms)')
                ax.set_ylabel('Density')
                fig.savefig('{}/ISI_analysis/ISImeanhistplot_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('mean ISI (ms)')
                ax.set_ylabel('Density')
                ax.set_xlim(0, np.mean(ISI_mean))
                fig.savefig('{}/ISI_analysis/ISImeanhistplot0_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('mean ISI (ms)')
                ax.set_ylabel('Density')
                ax.set_xlim(0, np.mean(ISI_mean) + 1*np.std(ISI_mean))
                fig.savefig('{}/ISI_analysis/ISImeanhistplot1_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('mean ISI (ms)')
                ax.set_ylabel('Density')
                ax.set_xlim(0, np.mean(ISI_mean) + 2*np.std(ISI_mean))
                fig.savefig('{}/ISI_analysis/ISImeanhistplot2_{}.png'.format(sim_dir, curr_idc))
            
            
            else:
                with open(ISI_file, 'a') as f:
                    print('Analysis:', analysand,end='\n\n', file=f)
                    
                    print('Spiking neurons (n spikes >= {}): no neurons'.format(min_sp_num), end='\n\n', file=f)
                    
            print('REPORT: ISI analysis concluded\n')
            
            
        return spiketime_list, zerolagCC_list, ISImean_list, ISIstd_list, ISICV_list, binned_spiketimes_list, return_expecvalue_list, return_var_list, autocorrvalue_list, autocorrt_list
       
    def DAcorrelation_analysis(self, analysis, sim_dir):
        
        print('REPORT: Starting DA correlation analysis\n')
        if not os.path.isdir('{}/correlationDA_analysis'.format(sim_dir)):
            os.mkdir('{}/correlationDA_analysis'.format(sim_dir))
        
       
        zerolagCC_list = []
        autocorrvalue_list = []
        autocorrt_list = []
     
        for analysand in analysis:
            curr_idc = 0
            print('Analysis {}'.format(curr_idc))
            ISI_file = '{}/correlationDA_analysis/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
           
            while os.path.isfile(ISI_file):
                curr_idc+= 1
                ISI_file = '{}/correlationDA_analysis/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
                
            group_n_stripes = analysand['Group']
            start_time = analysand['Start']
            stop_time = analysand['Stop']
            time_bin = analysand['Time bin']
            min_sp_num = analysand['Minimum spike number']
                
            neurons = self.group_name_to_neuron_idc(group_n_stripes)
            idc_duration = np.where((self.spikemonitor.t/ms>= start_time)&(self.spikemonitor.t/ms<stop_time))[0]
            neurons_duration, counts = np.unique(self.spikemonitor.i[idc_duration], return_counts=True)
            
            spiking = neurons_duration[np.where(counts>= min_sp_num)]
            neurons_idc = np.intersect1d(neurons, spiking)
            
            if len(neurons_idc) > 0:
                time_interval = stop_time - start_time
                time_interval_in_s = time_interval/1000
                time_bin_in_s = time_bin/1000
                N_bins = int(round(time_interval/time_bin, 0))
                
                expecvalue_list = []
                var_list = []
                r_list = []
                
     
                CC_arr = []
                autoC_arr = np.zeros((len(neurons_idc), N_bins))
                autoCt_arr = np.arange(N_bins)*time_bin
                
                totalneurons = len(neurons_idc)*(len(neurons_idc) - 1)//2
                k=0
                dez = 0
                for i in range(len(neurons_idc)):
                    for j in range(i):
                        k+= 1
                        perc = k/totalneurons*100
                        if perc//10 > dez:
                            dez = perc//10
                            print('{}% of cross-correlation concluded'.format(int(round(perc, 0))))
                        
                        T_spike1 = self.spikemonitor.t[self.spikemonitor.i==neurons_idc[i]]/ms
                        T_spike2 = self.spikemonitor.t[self.spikemonitor.i==neurons_idc[j]]/ms
                    
                        diff_spiketimes = []
                        for T in T_spike1:
                            diff_spiketimes.extend(T_spike2 - T)
                    
                        diff_spiketimes = np.asarray(diff_spiketimes)
                        
                        N0 = len(np.where((diff_spiketimes>-time_bin/2) & (diff_spiketimes<time_bin/2))[0])
                    
                        CC_arr.append(N0/time_interval_in_s - len(T_spike1) * len(T_spike2) * time_bin_in_s/time_interval_in_s**2)
                
                print('Cross-correlation concluded\n')
                
                k = 0
                dez = 0
                totalneurons = len(neurons_idc)
                for i in range(len(neurons_idc)):
                    perc = k/totalneurons*100
                    k+=1
                    if perc//10 > dez:
                        dez = perc//10
                        print('{}% of autocorrelation concluded'.format(int(round(perc, 0))))
                    T_spike = self.spikemonitor.t[self.spikemonitor.i==neurons_idc[i]]/ms
                    
                    diff_spiketimes = []
                    for T in T_spike:
                        diff_spiketimes.extend(T_spike - T)
                
                    diff_spiketimes = np.asarray(diff_spiketimes)
                    
                    diff_m = np.floor(diff_spiketimes/time_bin +1/2)
    
                    m_value, m_count = np.unique(diff_m, return_counts=True)
                    m_count = m_count[np.where((m_value>=0)&(m_value<N_bins))]
                    m_value = m_value[np.where((m_value>=0)&(m_value<N_bins))].astype(int)
                    
                    autoC_arr[i, m_value] = m_count
                    autoC_arr[i, :] = autoC_arr[i,:]/time_interval_in_s - len(T_spike)**2*time_bin_in_s/time_interval_in_s**2
                
                
                
                print('Autocorrelation concluded')
                
                with open(ISI_file, 'a') as f:
                    print('Analysis:', analysand,end='\n\n', file=f)
                    
                    print('Spiking neurons (n spikes >= {}):'.format(min_sp_num), len(neurons_idc), end='\n\n', file=f)
                      
                    print('CC(0) mean:', np.mean(CC_arr), file=f)
                    print('CC(0) std:', np.std(CC_arr), end='\n\n', file=f)
                 
                zerolagCC_list.append(CC_arr)
                autocorrvalue_list.append(autoC_arr)
                autocorrt_list.append(autoCt_arr)
                
                fig, ax = subplots()
                meanautoC = np.mean(autoC_arr, axis=0)
                ax.set_xlabel('Lag (ms)')
                ax.set_ylabel('Autocorrelation (Hz)')
                ax.plot(autoCt_arr[:len(autoCt_arr//4)], meanautoC[:len(autoCt_arr//4)])
                fig.savefig('{}/correlationDA_analysis/meanautoCtrace_{}.png'.format(sim_dir, curr_idc))
                
                autoCt_arr = np.asarray(autoCt_arr) 
                t100ms = autoCt_arr[autoCt_arr<300]
                fig, ax = subplots()
                ax.plot(autoCt_arr[:len(t100ms)], meanautoC[:len(t100ms)])
                ax.set_xlabel('Lag (ms)')
                ax.set_ylabel('Autocorrelation (Hz)')
                fig.savefig('{}/correlationDA_analysis/meanautoCtrace300ms_{}.png'.format(sim_dir, curr_idc))
    
                autoCt_arr = np.asarray(autoCt_arr) 
                t10_100ms = np.where((autoCt_arr>=2) & (autoCt_arr<300))[0]
                fig, ax = subplots()
                ax.plot(autoCt_arr[t10_100ms], meanautoC[t10_100ms])
                ax.set_xlabel('Lag (ms)')
                ax.set_ylabel('Autocorrelation (Hz)')
                fig.savefig('{}/correlationDA_analysis/meanautoCtrace2_300ms_{}.png'.format(sim_dir, curr_idc))
                
                autoCt_arr = np.asarray(autoCt_arr) 
                t10_500ms = np.where((autoCt_arr>=2) & (autoCt_arr<1000))[0]
                fig, ax = subplots()
                ax.plot(autoCt_arr[t10_500ms], meanautoC[t10_500ms])
                ax.set_xlabel('Lag (ms)')
                ax.set_ylabel('Autocorrelation (Hz)')
                fig.savefig('{}/correlationDA_analysis/meanautoCtrace2_1000ms_{}.png'.format(sim_dir, curr_idc))
                
                bins_lims, dens_lims, dens_arr = build_hist(CC_arr)
                
                fig, ax = subplots()
                ax.hist(CC_arr, bins=int(round(np.sqrt(len(CC_arr)), 0)))
                ax.set_xlabel('Zero-lag cross-correlation (Hz)')
                ax.set_ylabel('Number of pairs')
                fig.savefig('{}/correlationDA_analysis/crosscorhist_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.hist(CC_arr, bins=bins_lims, density=True)
                ax.set_xlabel('Zero-lag cross-correlation (Hz)')
                ax.set_ylabel('Density')
                fig.savefig('{}/correlationDA_analysis/crosscorhistdens_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Zero-lag cross-correlation (Hz)')
                ax.set_ylabel('Density')
                fig.savefig('{}/correlationDA_analysis/crosscorhistplot_{}.png'.format(sim_dir, curr_idc))
                
    
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Zero-lag cross-correlation (Hz)')
                ax.set_ylabel('Density')
                ax.set_xlim(np.mean(CC_arr) - 1* np.std(CC_arr), np.mean(CC_arr) + 1*np.std(CC_arr))
                fig.savefig('{}/correlationDA_analysis/crosscorhistplot1_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                ax.plot(dens_lims, dens_arr)
                ax.set_xlabel('Zero-lag cross-correlation (Hz)')
                ax.set_ylabel('Density')
                ax.set_xlim(np.mean(CC_arr) - 2* np.std(CC_arr), np.mean(CC_arr) + 2*np.std(CC_arr))
                fig.savefig('{}/correlationDA_analysis/crosscorhistplot2_{}.png'.format(sim_dir, curr_idc))
                
            
                
            else:
                with open(ISI_file, 'a') as f:
                    print('Analysis:', analysand,end='\n\n', file=f)
                    
                    print('Spiking neurons (n spikes >= {}): no neurons'.format(min_sp_num), end='\n\n', file=f)
                      
                
        print('REPORT: DA correlation analysis concluded\n')
            
        return zerolagCC_list, autocorrvalue_list, autocorrt_list
       
        
    def V_analysis(self, analysis, sim_dir):
            
        print('REPORT: Starting V analysis\n')
        
        if 'V' not in self.neuronmonitors:
            print("ERROR: No 'V' monitor\n")
            return
        
        if not os.path.isdir('{}/V_analysis'.format(sim_dir)):
            os.mkdir('{}/V_analysis'.format(sim_dir))
                       
        Vmean_list = []
        Vstd_list = []
        Vsubthres_list = []
               
        monitor = self.neuronmonitors['V']
        monitor_t = monitor.t/ms
        
        for analysand in analysis:
            curr_idc = 0
            V_file = '{}/V_analysis/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
        
            while os.path.isfile(V_file):
                curr_idc+= 1
                V_file = '{}/V_analysis/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
            
            group = analysand['Group']
            start_time = analysand['Start']
            stop_time = analysand['Stop']
            min_sp = analysand['Minimum spike number']
            
            neuron_idc = self.group_name_to_neuron_idc(group)
                    
            if len(np.setdiff1d(neuron_idc, monitor.record))>0:
                print('WARNING: Skipping current group because there are neurons not recorded')
                break    
        
            spiking = np.where(self.spikemonitor.count[neuron_idc] >= min_sp)[0] 
            
            V_up_all = self.NeuPar[4, spiking]
            V_r_all = self.NeuPar[8, spiking]
            V_T_all = self.NeuPar[9, spiking]
            
            N = len(spiking)
            monitor_t = np.where((monitor.t/ms >= start_time) & (monitor.t/ms < stop_time))[0]
            V_array = monitor.V[spiking, :]/mV
            V_trimmed = V_array[:, monitor_t]
            
            Vindmean = []
            Vindstd = []
            Vsubthres = []
            
            for n in range(N):
                V1 = V_trimmed[n, :]
                V_up = V_up_all[n]
                V_r = V_r_all[n]
                V_T = V_T_all[n]
                
                Vex = extract_spikes(V1, V_up, V_r, V_T, timestep=self.time_step)
                Vexmean = np.mean(Vex)
                Vindmean.append(Vexmean)
                Vindstd.append(np.std(Vex))
                Vsubthres.append(V_T - Vexmean)
                
             
            with open(V_file, 'a') as f:
                print('Analysis:', analysand,end='\n\n', file=f)
                
                print('Spiking neurons (n spikes >= {}):'.format(min_sp), N, end='\n\n', file=f)
                
                print('Individual Vstd mean:', np.mean(Vindstd), file=f)
                print('Individual Vstd std:', np.std(Vindstd), end='\n', file=f)
                print('Individual Vstd min:', np.min(Vindstd), end='\n', file=f)
                print('Individual Vstd max:', np.max(Vindstd), end='\n\n', file=f)
                
                print('Individual Vmean mean:', np.mean(Vindmean), file=f)
                print('Individual Vmean std:', np.std(Vindmean), end='\n', file=f)
                print('Individual Vmean min:', np.min(Vindmean), end='\n', file=f)
                print('Individual Vmean max:', np.max(Vindmean), end='\n\n', file=f)
                
                print('Individual V_T - Vmean mean:', np.mean(Vsubthres), file=f)
                print('Individual V_T - Vmean std:', np.std(Vsubthres), end='\n', file=f)
                print('Individual V_T - Vmean min:', np.min(Vsubthres), end='\n', file=f)
                print('Individual V_T - Vmean max:', np.max(Vsubthres), end='\n\n', file=f)
                                
                  
            bins_lims, dens_lims, dens_arr = build_hist(Vindstd)
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Number of neurons')
            ax.hist(Vindstd, bins=int(round(np.sqrt(len(Vindstd)), 0)))
            fig.savefig('{}/V_analysis/Vindstd_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Density')
            ax.hist(Vindstd, bins=bins_lims, density=True)
            fig.savefig('{}/V_analysis/Vindstddens_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            fig.savefig('{}/V_analysis/Vindstdplot_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0, np.mean(Vindstd))
            fig.savefig('{}/V_analysis/Vindstdplot0_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0, np.mean(Vindstd) + 1 * np.std(Vindstd))
            fig.savefig('{}/V_analysis/Vindstdplot1_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0,  np.mean(Vindstd) + 2* np.std(Vindstd))
            fig.savefig('{}/V_analysis/Vindstdplot2_{}.png'.format(sim_dir, curr_idc))
            
            
            
            bins_lims, dens_lims, dens_arr = build_hist(Vindmean)
            
            fig, ax = subplots()
            ax.set_xlabel('V mean (mV)')
            ax.set_ylabel('Number of neurons')
            ax.hist(Vindmean, bins=int(round(np.sqrt(len(Vindmean)), 0)))
            fig.savefig('{}/V_analysis/Vindmean_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V mean (mV)')
            ax.set_ylabel('Density')
            ax.hist(Vindmean, bins=bins_lims, density=True)
            fig.savefig('{}/V_analysis/Vindmeandens_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V mean (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            fig.savefig('{}/V_analysis/Vindmeanplot_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V mean (mV)')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0, np.mean(Vindmean))
            fig.savefig('{}/V_analysis/Vindmeanplot0_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V mean (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0, np.mean(Vindmean) + 1 * np.std(Vindmean))
            fig.savefig('{}/V_analysis/Vindmeanplot1_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V mean (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0,  np.mean(Vindmean) + 2* np.std(Vindmean))
            fig.savefig('{}/V_analysis/Vindmeanplot2_{}.png'.format(sim_dir, curr_idc))
                       
            
            bins_lims, dens_lims, dens_arr = build_hist(Vsubthres)
            
            fig, ax = subplots()
            ax.set_xlabel('V_T - V mean (mV)')
            ax.set_ylabel('Number of neurons')
            ax.hist(Vsubthres, bins=int(round(np.sqrt(len(Vsubthres)), 0)))
            fig.savefig('{}/V_analysis/Vsubthres_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V_T - V mean (mV)')
            ax.set_ylabel('Density')
            ax.hist(Vsubthres, bins=bins_lims, density=True)
            fig.savefig('{}/V_analysis/Vsubthresdens_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V_T - V mean (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            fig.savefig('{}/V_analysis/Vsubthresplot_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V_T - V mean (mV)')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0, np.mean(Vsubthres))
            fig.savefig('{}/V_analysis/Vsubthresplot0_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V_T - V mean (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0, np.mean(Vsubthres) + 1 * np.std(Vsubthres))
            fig.savefig('{}/V_analysis/Vsubthresplot1_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V_T - V mean (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0,  np.mean(Vsubthres) + 2* np.std(Vsubthres))
            fig.savefig('{}/V_analysis/Vsubthresplot2_{}.png'.format(sim_dir, curr_idc))
            
            Vmean_list.append(Vindmean)
            Vstd_list.append(Vindstd)
            Vsubthres_list.append(Vsubthres)
            
        print('REPORT:V analysis concluded\n')
        return Vmean_list, Vstd_list, Vsubthres_list

    def Vcorr_analysis(self, analysis, sim_dir):
            
        print('REPORT: Starting Vcorr analysis\n')
        
        if 'V' not in self.neuronmonitors:
            print("ERROR: No 'V' monitor\n")
            return
        
        if not os.path.isdir('{}/Vcorr_analysis'.format(sim_dir)):
            os.mkdir('{}/Vcorr_analysis'.format(sim_dir))
                       
        Vgroup_list = []
        VzerolagCC_list = []
        VnormalizedzerolagCC_list = []
        Vindividualstd_list = []
        Vindividualmean_list = []
               
        monitor = self.neuronmonitors['V']
        monitor_t = monitor.t/ms
        group_list = []
        
        for analysand in analysis:
            curr_idc = 0
            V_file = '{}/Vcorr_analysis/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
        
            while os.path.isfile(V_file):
                curr_idc+= 1
                V_file = '{}/Vcorr_analysis/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
            
            group = analysand['Group']
            start_time = analysand['Start']
            stop_time = analysand['Stop']
            min_sp = analysand['Minimum spike number']
            
            neuron_idc = self.group_name_to_neuron_idc(group)
                    
            if len(np.setdiff1d(neuron_idc, monitor.record))>0:
                print('WARNING: Skipping current group because there are neurons not recorded')
                break    
        
            spiking = np.where(self.spikemonitor.count[neuron_idc] >= min_sp)[0] 
            duration = stop_time-start_time
            
            N = len(spiking)
            monitor_t = np.where((monitor.t/ms >= start_time) & (monitor.t/ms <= stop_time))[0]
            V_array = monitor.V[spiking, :]/mV
            V_trimmed = V_array[:, monitor_t]
            
            V_group = np.sum(V_trimmed, axis=0)/N
            sig_pop = np.std(V_group)**2
            
            V_var = np.zeros(N)
            
            for i in range(N):
                V_var[i] = np.std(V_trimmed[i, :])**2    
                
            qui = sig_pop/((1/N)*np.sum(V_var))
            
            Vm_arr = np.zeros((N, len(monitor_t)))
            Vstd_arr = np.zeros(N)
            Vmean_arr = np.zeros(N)
            
            for n in range(N):
                V1 = V_trimmed[n, :]
                             
                V1m = np.mean(V1)
                Vm_arr[n, :] = V1-V1m
                Vstd_arr[n] = np.std(V1)
                Vmean_arr[n] = V1m
                
            V_corr = []
            
            for i in range(N):
                for j in range(i):
                    _ = np.correlate(Vm_arr[i, :], Vm_arr[j, :])[0]*self.time_step/duration
                    V_corr.append(_)
                    
            V_norm_corr = []
            k=0
            
            for i in range(N):
                for j in range(i):
                    _ = V_corr[k]/(Vstd_arr[i] * Vstd_arr[j])
                    V_norm_corr.append(_)
                    k+=1        
              
            with open(V_file, 'a') as f:
                print('Analysis:', analysand,end='\n\n', file=f)
                
                print('Spiking neurons (n spikes >= {}):'.format(min_sp), N, end='\n\n', file=f)
                
                print('Individual Vstd mean:', np.mean(Vstd_arr), file=f)
                print('Individual Vstd std:', np.std(Vstd_arr), end='\n', file=f)
                print('Individual Vstd min:', np.min(Vstd_arr), end='\n', file=f)
                print('Individual Vstd max:', np.max(Vstd_arr), end='\n\n', file=f)
                
                
                
                print('Mean group V:', np.mean(V_group), file=f)
                print('Std group V:', sig_pop, file=f)
                print('Min group V:', np.min(V_group), file=f)
                print('Max group V:', np.max(V_group), file=f)
                print('Qui:', qui, end='\n\n', file=f)
                
                print('Zero-lag cross-correlation', file=f)
                print('Mean:', np.mean(V_corr), file=f)
                print('Std:', np.std(V_corr), end='\n\n', file=f)
                
                print('Normalized zero-lag cross-correlation', file=f)
                print('Mean:', np.mean(V_norm_corr), file=f)
                print('Std:', np.std(V_norm_corr), end='\n\n', file=f)
                print('-'*50, end='\n\n', file=f)

                Vgroup_list.append(V_group)
                VzerolagCC_list.append(V_corr)
                VnormalizedzerolagCC_list.append(V_norm_corr)
                Vindividualstd_list.append(Vstd_arr)
                Vindividualmean_list.append(Vmean_arr)
                
            fig, ax = subplots()
            ax.plot(monitor_t, V_group)
            ax.set_xlabel('t (ms)')
            ax.set_ylabel('V (mV)')
        
            fig.savefig('{}/Vcorr_analysis/V_group_{}.png'.format(sim_dir, curr_idc))
            
            bins_lims, dens_lims, dens_arr = build_hist(V_corr)
            
            fig, ax = subplots()
            ax.hist(V_corr, bins=int(round(np.sqrt(len(V_corr)), 0)))
            ax.set_xlabel('Zero-lag cross-correlation (mV^2)')
            ax.set_ylabel('Number of pairs')
            fig.savefig('{}/Vcorr_analysis/V_curr_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.hist(V_corr, bins=bins_lims, density=True)
            ax.set_xlabel('Zero-lag cross-correlation (mV^2)')
            ax.set_ylabel('Density')
            fig.savefig('{}/Vcorr_analysis/V_currdens_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.plot(dens_lims, dens_arr)
            ax.set_xlabel('Zero-lag cross-correlation (mV^2)')
            ax.set_ylabel('Density')
            fig.savefig('{}/Vcorr_analysis/V_currplot_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.plot(dens_lims, dens_arr)
            ax.set_xlabel('Zero-lag cross-correlation (mV^2)')
            ax.set_ylabel('Density')
            ax.set_xlim(np.mean(V_corr) - 1* np.std(V_corr), np.mean(V_corr) + 1* np.std(V_corr))
            fig.savefig('{}/Vcorr_analysis/V_currplot1_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.plot(dens_lims, dens_arr)
            ax.set_xlabel('Zero-lag cross-correlation (mV^2)')
            ax.set_ylabel('Density')
            ax.set_xlim(np.mean(V_corr) - 2* np.std(V_corr), np.mean(V_corr) + 2* np.std(V_corr))
            fig.savefig('{}/Vcorr_analysis/V_currplot2_{}.png'.format(sim_dir, curr_idc))
            
            
            
            bins_lims, dens_lims, dens_arr = build_hist(V_norm_corr)
            
            fig, ax = subplots()
            ax.set_xlabel('Zero-lag cross-correlation')
            ax.set_ylabel('Number of pairs')
            ax.hist(V_norm_corr, bins=int(round(np.sqrt(len(V_norm_corr)), 0)))
            fig.savefig('{}/Vcorr_analysis/V_norm_curr_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('Zero-lag cross-correlation')
            ax.set_ylabel('Density')
            ax.hist(V_norm_corr, bins=bins_lims, density=True)
            fig.savefig('{}/Vcorr_analysis/V_norm_currdens_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('Zero-lag cross-correlation')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            fig.savefig('{}/Vcorr_analysis/V_norm_currplot_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('Zero-lag cross-correlation')
            ax.set_ylabel('Density')
            ax.set_xlim(np.mean(V_norm_corr) - 1*np.std(V_norm_corr), np.mean(V_norm_corr) - 1*np.std(V_norm_corr))
            ax.plot(dens_lims, dens_arr)
            fig.savefig('{}/Vcorr_analysis/V_norm_currplot1_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('Zero-lag cross-correlation')
            ax.set_ylabel('Density')
            ax.set_xlim(np.mean(V_norm_corr) - 2*np.std(V_norm_corr), np.mean(V_norm_corr) - 2*np.std(V_norm_corr))
            ax.plot(dens_lims, dens_arr)
            fig.savefig('{}/Vcorr_analysis/V_norm_currplot2_{}.png'.format(sim_dir, curr_idc))
            
            bins_lims, dens_lims, dens_arr = build_hist(Vstd_arr)
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Number of neurons')
            ax.hist(Vstd_arr, bins=int(round(np.sqrt(len(Vstd_arr)), 0)))
            fig.savefig('{}/Vcorr_analysis/Vstd_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Density')
            ax.hist(Vstd_arr, bins=bins_lims, density=True)
            fig.savefig('{}/Vcorr_analysis/Vstddens_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            fig.savefig('{}/Vcorr_analysis/Vstdplot_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0, np.mean(Vstd_arr))
            fig.savefig('{}/Vcorr_analysis/Vstdplot0_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0, np.mean(Vstd_arr) + 1 * np.std(Vstd_arr))
            fig.savefig('{}/Vcorr_analysis/Vstdplot1_{}.png'.format(sim_dir, curr_idc))
            
            fig, ax = subplots()
            ax.set_xlabel('V standard deviation (mV)')
            ax.set_ylabel('Density')
            ax.plot(dens_lims, dens_arr)
            ax.set_xlim(0,  np.mean(Vstd_arr) + 2* np.std(Vstd_arr))
            fig.savefig('{}/Vcorr_analysis/Vstdplot2_{}.png'.format(sim_dir, curr_idc))
                       
            
        print('REPORT:Vcorr analysis concluded\n')
        return monitor_t,Vindividualstd_list, Vindividualmean_list, Vgroup_list, VzerolagCC_list, VnormalizedzerolagCC_list

    
               
    def frequence_spectrum(self, analysis, sim_dir):

        print('REPORT: Starting frequence analysis\n')
        
        Imonitort_list = []
        LFPfrequence_list = []
        I_list = []
        LFPpower_list = []
        MALFPfrequence_list = []
        MALFPpower_list = []
        
        for analysand in analysis:
            
            source = analysand['Source']
            
            if source not in self.neuronmonitors:
                print("ERROR: No '{}' monitor\n".format(source))
                break
            
            if not os.path.isdir('{}/{}_frequence'.format(sim_dir, source)):
                os.mkdir('{}/{}_frequence'.format(sim_dir, source))
                
            monitor = self.neuronmonitors[source]
            curr_idc = 0
            Fq_file = '{}/{}_frequence/Report_{}_{}.txt'.format(sim_dir, source,curr_idc, sim_dir)
            
            while os.path.isfile(Fq_file):
                curr_idc+= 1
                Fq_file = '{}/{}_frequence/Report_{}_{}.txt'.format(sim_dir, source, curr_idc, sim_dir)
            
            group = analysand['Group']
            start_time = analysand['Start']
            stop_time = analysand['Stop']
            min_sp = analysand['Minimum spike number']
            max_sp = analysand['Maximum spike number']
            
            neuron_idc = self.group_name_to_neuron_idc(group)
            count_idc = np.where((self.spikemonitor.count>= min_sp)&(self.spikemonitor.count<max_sp))[0]
            
            neuron_idc = np.intersect1d(neuron_idc, count_idc)
        
            I_arr = eval('monitor.{}[neuron_idc]/pA'.format(source))
            
            t_arr = monitor.t/ms
            idc_arr = np.where((t_arr>=start_time) & (t_arr<stop_time))[0]
            
            I_arr = I_arr[:, idc_arr]
            t_arr = t_arr[idc_arr]

            I_arr = np.sum(I_arr, axis=0)
            fs = 1000/self.time_step
            
            freq_arr, power_arr = periodogram(I_arr, fs)
            
            Imonitort_list.append(t_arr)
            LFPfrequence_list.append(freq_arr)
            I_list.append(I_arr)
            LFPpower_list.append(power_arr)
            
            log_power = np.log(power_arr)
            log_freq = np.log(freq_arr)
            
            ma_log_power = moving_average(log_power, 11)
            ma_log_freq = log_freq[5:-5]
            
            MALFPfrequence_list.append(ma_log_power)
            MALFPpower_list.append(ma_log_freq)
            
            with open(Fq_file, 'a') as f:
                print()
            
                print('Analysis:', analysand,end='\n\n', file=f)
                
                print('Spiking neurons ({} <= n spikes <{} ):'.format(min_sp, max_sp), len(neuron_idc), end='\n\n', file=f)
            
                print('Population {}\n', file=f)
                print('Min pop {}:'.format(source), np.min(I_arr), file=f)
                print('Max pop {}:'.format(source), np.max(I_arr), file=f)
                print('Mean pop {}:'.format(source), np.mean(I_arr), file=f)
                print('Std pop {}:'.format(source), np.std(I_arr), end='\n\n', file=f)
            
                print('Mean individual {}:'.format(source), np.mean(I_arr)/len(count_idc), file=f)    
            
            fig, ax = subplots()
            xlim(1, 7)

            ax.plot(log_freq, log_power)
            ax.set_xlabel('log(Frequency) (log[Hz])')
            ax.set_ylabel('log(Power) (arbitrary unit)')
            plt.savefig('{}/{}_frequence/Fqspectrum_{}_{}.png'.format(sim_dir, source, source, curr_idc))
            
            fig, ax = subplots()
            xlim(1, 7)
            ax.set_xlabel('log(Frequency) (log[Hz])')
            ax.set_ylabel('log(Power) (arbitrary unit)')
            ax.plot(ma_log_freq, ma_log_power)
            plt.savefig('{}/{}_frequence/Fqspectrum_ma_{}_{}.png'.format(sim_dir,source, source, curr_idc))
            
            fig, ax = subplots()
            ax.plot(t_arr, I_arr)
            ax.set_xlabel('t (ms)')
            ax.set_ylabel('I_{} (pA)'.format(source))
            plt.savefig('{}/{}_frequence/{}_{}.png'.format(sim_dir, source, source, curr_idc))           
            
        print('REPORT: Frequency analysis concluded\n')
        
        return Imonitort_list, I_list, LFPfrequence_list, LFPpower_list, MALFPfrequence_list, MALFPpower_list
        
    def population_rate(self, analysis, sim_dir):
        
        if not os.path.isdir('{}/Rate'.format(sim_dir)):
            os.mkdir('{}/Rate'.format(sim_dir))
                
        popratet_lists = []
        popratecount_lists = []
        popratefreq_lists = []
        
        
        for analysand in analysis:
            
            curr_idc = 0
            Rate_file = '{}/Rate/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
            while os.path.isfile(Rate_file):
                curr_idc+= 1
                Rate_file = '{}/Rate/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
                     
            group = analysand['Group']
            start = analysand['Start']
            stop = analysand['Stop']
            time_bins = analysand['Time bin']
            MA = analysand['Moving average']
            duration = (stop-start)/1000
            
            neuron_idc = self.group_name_to_neuron_idc(group)
            
            if len(neuron_idc) > 0:
                idc_arr = np.isin(self.spikemonitor.i, neuron_idc)
                
                t_arr = self.spikemonitor.t[idc_arr]/ms
                
                t_arr = t_arr[np.where((t_arr>=start)&(t_arr<stop))[0]]
                
                N_bins = int(round((stop - start)/time_bins, 0))
                count_arr = []
                time_arr = []
                fq_arr = []
                
                for i in range(N_bins):
                    t0 = start + i*time_bins
                    t1 = start + (i+1)*time_bins
                    time_arr.append((t1+t0)/2)
                    count = len(np.where((t_arr>=t0)&(t_arr<t1))[0])
                    count_arr.append(count)
                    fq_arr.append((count/len(neuron_idc))*(1000/time_bins))
                    
                popratet_lists.append(time_arr)
                popratecount_lists.append(count_arr)
                popratefreq_lists.append(fq_arr)
                
                with open(Rate_file, 'a') as f:
                
                    print('Analysis:', analysand,end='\n\n', file=f)
                    
                    print('Number of neurons: {}'.format(len(neuron_idc)), file=f)
                    print('Total count: {}'.format(np.sum(count_arr)),file=f)
                    print('Mean spikes per neurons: {}\n'.format(round(np.sum(count_arr)/len(neuron_idc),2)), file=f)
                    
                    print('Population overall frequence: {} Hz'.format(round(np.sum(count_arr)/(len(neuron_idc)*duration),2)), file=f)
                    
                    print('Population frequence  mean: {} Hz'.format(round(np.mean(fq_arr),2)), file=f)
                    print('Population frequence std: {} Hz\n'.format(round(np.std(fq_arr),2)), file=f)
                    
                    for i in range(len(time_arr)):
                        print('Time: {} ms'.format(time_arr[i]), file=f)
                        print('Count: {}'.format(count_arr[i]), file=f)
                        print('Mean frequence: {}\n'.format(round(fq_arr[i],2)), file=f)
                                   
                fig, ax = subplots()
                ax.plot(time_arr, fq_arr)
                ax.set_xlabel('time (ms)')
                ax.set_ylabel('Firing rate (Hz)')
                fig.savefig('{}/Rate/Pop_rate_{}.png'.format(sim_dir, curr_idc))
                
                fig, ax = subplots()
                
                if MA >1:
                    MAfq_arr = moving_average(fq_arr, MA)
                    lim = (len(fq_arr) - len(MAfq_arr))//2
                    ax.plot(time_arr[lim:-lim], MAfq_arr)
                    ax.set_xlabel('time (ms)')
                    ax.set_ylabel('Firing rate (Hz)')
                    fig.savefig('{}/Rate/MAPop_rate_{}.png'.format(sim_dir, curr_idc))
                    
            else:
                with open(Rate_file, 'a') as f:
                
                    print('Analysis:', analysand,end='\n\n', file=f)
                    
                    print('Number of neurons: 0', file=f)
                
        return popratet_lists, popratecount_lists, popratefreq_lists
        
    def rate_distribution(self, analysis, sim_dir):
        
        ratedistribution_total_list = []
        ratedistribution_count_list = []
        ratedistribution_neuron_list = []
        
        for analysand in analysis:
                   
            if not os.path.isdir('{}/RateDistribution'.format(sim_dir)):
                os.mkdir('{}/RateDistribution'.format(sim_dir))

            curr_idc = 0
            Fq_file = '{}/RateDistribution/Report_{}_{}.txt'.format(sim_dir,curr_idc, sim_dir)
           
            while os.path.isfile(Fq_file):
                curr_idc+= 1
                Fq_file = '{}/RateDistribution/Report_{}_{}.txt'.format(sim_dir, curr_idc, sim_dir)
            
            group = analysand['Group']
            start_time = analysand['Start']
            stop_time = analysand['Stop']
            bins = analysand['Bins']    
            duration = (stop_time-start_time)/1000
            
            neuron_idc = self.group_name_to_neuron_idc(group)
            
            
            Nneurons = len(neuron_idc)
            
            
            count_list = []
            neuron_list = []
            
            if Nneurons:
                counts = self.spikemonitor.count[neuron_idc]
                fqs = counts/duration
                Nspikes = sum(counts)
                
                distr = []
                lower = 0
                bins.sort()
                
                ratedistribution_total_list.append(Nneurons)
                
                with open(Fq_file, 'a') as f:
                    
                    print("Group: {}".format(group), end='\n\n', file = f)
                    print('Number of neurons: {}'.format(Nneurons), file = f)
                    print('Mean frequency: {}'.format(np.mean(fqs)), end='\n\n', file=f)
                                   
                    for i in range(len(bins)):
                        if i == 0:
                            Nfq = len(np.where((fqs>= lower) & (fqs <= bins[i]))[0])
                            NeuronsFq = neuron_idc[np.where((fqs>= lower) & (fqs <= bins[i]))[0]]
                            count_list.append(Nfq)
                            neuron_list.append(NeuronsFq)
                            text = "0 Hz <= f <= {} Hz: {}% ({} neurons)".format(round(bins[i], 2), round(100*Nfq/Nneurons, 2), Nfq)
                            print(text, file=f)
                            lower = bins[i]
                        else:
                            Nfq = len(np.where((fqs> lower) & (fqs <= bins[i]))[0])
                            NeuronsFq = neuron_idc[np.where((fqs> lower) & (fqs <= bins[i]))[0]]
                            count_list.append(Nfq)
                            neuron_list.append(NeuronsFq)
                            text = "{} Hz < f <= {} Hz: {}% ({} neurons)".format(round(lower, 2), round(bins[i], 2), round(100*Nfq/Nneurons, 2), Nfq)
                            print(text, file=f)
                            lower = bins[i]
                    else:
                        Nfq = len(np.where((fqs > lower))[0])
                        NeuronsFq = neuron_idc[np.where((fqs > lower))[0]]
                        count_list.append(Nfq)
                        neuron_list.append(NeuronsFq)
                        text="{}Hz < f: {}% ({} neurons)".format(round(lower,2), round(100*Nfq/Nneurons, 2), Nfq)
                        print(text, file=f, end='\n'+'-'*40)
            
            else:
                with open(Fq_file, 'a') as f:
                    
                    print("Group: {}".format(group), end='\n\n', file = f)
                    print('Number of neurons: {}'.format(Nneurons), file = f)
                 
            ratedistribution_count_list.append(count_list)
            ratedistribution_neuron_list.append(neuron_list)
            
        return ratedistribution_total_list, ratedistribution_count_list, ratedistribution_neuron_list
                    
    def raster_plot(self, sim_dir, tlims=False):
        
        fig, ax = subplots()  
        
        if type(tlims) is not bool:
            xlim(tlims[0], tlims[1])
            
        if not os.path.isdir('{}'.format(sim_dir)):
            os.mkdir('{}'.format(sim_dir))
            
        ax.scatter(self.spikemonitor.t/ms, self.spikemonitor.i, 1)
        ax.set_xlabel('time (ms)')
        ax.set_ylabel('Neuron index')
        
        curr_idc = 0
        RPfile = '{}/Raster_plot_{}.png'.format(sim_dir, curr_idc)
        while os.path.isfile(RPfile):
            curr_idc += 1
            RPfile = '{}/Raster_plot_{}.png'.format(sim_dir, curr_idc)
        
        fig.savefig(RPfile)

    
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

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def build_hist(list_):
    
    arr_= np.asarray(list_)
    ntotal = np.sum(list_)
    nbinlimits = int(round(np.sqrt(len(arr_))))
    
    bins_list = []
    for i in np.linspace(0, 100, nbinlimits+1):
        bins_list.append(np.percentile(arr_, i))
    bins_arr = np.asarray(bins_list)
    
    bins_diff = np.diff(bins_arr)
    bins_mid = bins_arr[:-1] + bins_diff/2
    
    dens_arr, _ = np.histogram(arr_, bins=bins_arr, density=True)
    
    dens_arr = np.asarray([0,] + list(dens_arr) + [0,])
    bins_mid_limits = np.asarray([bins_list[0],] + list(bins_mid) + [bins_list[-1],])
    
    return bins_arr, bins_mid_limits, dens_arr
    
    

    
    