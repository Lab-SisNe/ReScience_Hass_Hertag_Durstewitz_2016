import brian2 as br2
import numpy as np
from BasicsSetup import BaseClass
from NetworkSetup import network_setup, Network
from numpy.matlib import repmat
from dataclasses import dataclass
import time as tm
import os, shutil
from Auxiliar import time_report
from BasicsSetup import GROUPS_SETS as SETS, MEMB_PARAM_NAMES


br2.BrianLogger.suppress_name('resolution_conflict', 'device')

###HERE
class Cortex(BaseClass):
    
    @time_report('Cortex setup')
    def setup(Ncells, Nstripes, constant_stimuli, method, dt,
              transient=0, basics_scales=None, fluctuant_stimuli=None, seed=None):
        network = network_setup(Ncells, Nstripes, basics_scales, seed)
        
        return __class__(network, constant_stimuli, method, dt,
                         transient, fluctuant_stimuli, seed)
    
    @time_report('Cortex setup (with network loading)')
    def load(path, constant_stimuli, method, dt,
             transient=0, fluctuant_stimuli=None):
        
        network = Network.load(path)
        seed = network.seed
        return __class__(network, constant_stimuli,  method, dt,
                         transient, fluctuant_stimuli, seed)
    

    def __init__(self, network, constant_stimuli, method, dt, 
                  transient=0, fluctuant_stimuli=None,  seed=None):
        
        np.random.seed(seed)
        br2.seed(seed)
        self.network = network
        self.method = method
        self.dt = dt * br2.ms 
        
        if constant_stimuli is not None and isinstance(constant_stimuli[-1], (int,float)):
            constant_stimuli = [constant_stimuli]
        self.constant_stimuli = constant_stimuli
        
        if fluctuant_stimuli is not None and isinstance(fluctuant_stimuli[-1], (int,float)):
            fluctuant_stimuli = [fluctuant_stimuli]
        self.fluctuant_stimuli = fluctuant_stimuli
        
        self.transient = transient
        self.seed = seed
            
        self.net = br2.Network()      
        self.set_fluctuant_stimuli()
        
        membrane_events = self.get_membrane_events_dict()     
        membr_model = network.basics.equations.membr_model.format(self.fluctuant_str)
        
        self.neurons = br2.NeuronGroup(N=network.basics.struct.Ncells_total,
                                       model=membr_model,
                                       threshold=network.basics.equations.membr_threshold, 
                                       reset=network.basics.equations.membr_reset,
                                       events=membrane_events,
                                       method=method, 
                                       refractory=5*br2.ms,
                                       dt=self.dt)
       
        self.set_constant_stimuli()
        self.set_membrane_events()
        self.set_neuron_params()
        self.set_auxiliar_vars()
        self.set_initial_state()
           
        syn_dict_pathway = self.get_syn_pathway_dict()
        self.synapses = br2.Synapses(self.neurons, self.neurons, model=network.basics.equations.syn_model, on_pre=syn_dict_pathway,\
                                 method=method, dt=self.dt)
        
        if self.network.syn_pairs.shape[1]>0:
            self.synapses.connect(i=self.network.syn_pairs[1,:], j=self.network.syn_pairs[0,:])
        else:
            self.synapses.active = False
        
        self.set_synapses()            
        self.spikemonitor = br2.SpikeMonitor(self.neurons)
        if self.transient>0:
            self.spikemonitor.active=False
        self.net.add(self.neurons, self.synapses, self.spikemonitor,
                     *self.event_monitors.values())
                    
        self.external_stimuli=NetworkHolder()
        
        self.neuron_monitors=NetworkHolder()
        self.synapse_monitors=NetworkHolder()
        
        self.longrun_neuron_monitors=NetworkHolder()
        self.longrun_synapse_monitors=NetworkHolder()
        self.recorded=VarHolder()
        
        self.monitor_schedule={}
        self.longrun_monitor_control=[]
    
    def save(self, path):
        self.network.save(path)
    
    @time_report('Cortex simulation')
    def run(self, t, erase_longrun=True):
        
        print_out = "::   REPORT: Preparing simulation   ::"
        
        print('.'*len(print_out))
        print(print_out)
        print('.'*len(print_out), end='\n\n')
        
        print('Time step: {} ms'.format(self.dt/br2.ms))
        print('Integration method:', self.method)
        print('Seed:', self.seed, end='\n\n')
            
        ### fluctuant_array needs to be declared here if fluctuant stimuli are set
        ### fluctuant_array appears inside the model code of NeuronGroup in self.neurons
        if self.fluctuant_array is not None:
            fluctuant_array = br2.TimedArray(self.fluctuant_array, dt=self.dt)
            
        schedule_times = np.asarray(list(self.monitor_schedule.keys()))
        t0 = self.neurons.t/br2.ms
             
        longrun_times = []
        longrun_all_times = []
        for longrun_control in self.longrun_monitor_control:
            longrun_stop = longrun_control['stop'] if longrun_control['stop'] is not None and longrun_control['stop']<=t else t      
            times = np.arange(longrun_control['start'], longrun_stop, longrun_control['interval'])
            times = np.append(times, longrun_stop)
            longrun_times.append(times)
            longrun_all_times.extend(times)
    
        current_schedule_times = np.concatenate((schedule_times, longrun_all_times, [self.transient,]))
        current_schedule_times = current_schedule_times[np.where((current_schedule_times>=t0)&(current_schedule_times<t0+t))]
        current_schedule_times = np.unique(current_schedule_times)
        
        if t0 not in current_schedule_times:
            current_schedule_times = np.append(current_schedule_times, t0)
        if t0+t not in current_schedule_times:
            current_schedule_times = np.append(current_schedule_times, t0+t)
        
        current_schedule_times = np.sort(current_schedule_times)
        current_schedule_intervals = np.diff(current_schedule_times)
        
        print('Start time: {} ms'.format(t0))
        print('Stop time: {} ms'.format(t0+t))
        print('Total duration : {:.2f} ms (from {:.2f} to {:.2f} ms)'.format(t, t0, t0+t), end='\n\n')
        print('Simulation points:', ', '.join(['{} ms'.format(i) for i in current_schedule_times]))
        print(len(current_schedule_intervals), 'simulation segment(s) of duration:', ', '.join(np.round(current_schedule_intervals,2).astype(str)), 'ms', end='\n\n')
        
        if t0 in schedule_times:
            for monitor in self.monitor_schedule[t0]['start']:
                monitor.active=True
            for monitor in self.monitor_schedule[t0]['stop']:
                monitor.active=False
        
        for i in range(len(current_schedule_intervals)):
            print_out1a = '::   Starting segment {} out of {}'.format(i+1, len(current_schedule_intervals))
            print_out2a = '::   Duration {:.2f} ms (from {:.2f} to {:.2f} ms)'.format(current_schedule_intervals[i],current_schedule_times[i], current_schedule_times[i+1])
            print_out1 = print_out1a + ' '*(max(len(print_out1a), len(print_out2a))-len(print_out1a)) + '   ::'
            print_out2 = print_out2a + ' '*(max(len(print_out1a), len(print_out2a))-len(print_out2a)) + '   ::'

            print('.'*len(print_out1))
            print(print_out1)
            print(print_out2)
            print('.'*len(print_out1), end='\n\n')
            self.net.run(current_schedule_intervals[i]*br2.ms, report='text', report_period=10*br2.second)
            
            t1 = current_schedule_times[i+1]
            if t1 >= self.transient:
                self.spikemonitor.active=True
                
            if t1 in longrun_all_times or t1==current_schedule_times[-1]:

                for l in range(len(longrun_times)):
      
                    if t1 in longrun_times[l] or t1==current_schedule_times[-1]:
                        self.process_longrun(l, t1)
            
            if t1 in schedule_times:
                for monitor in self.monitor_schedule[t1]['start']:
                    monitor.active=True
                for monitor in self.monitor_schedule[t1]['stop']:
                    monitor.active=False
        print()
        self.restore_longrun(erase_longrun)
    
        
    def set_fluctuant_stimuli(self):
        if self.fluctuant_stimuli is not None:
            for target, function, start, stop in self.fluctuant_stimuli:
                neuron_idcs = self.neuron_idcs(target)
                
                Nsteps_total = round(stop/(self.dt/br2.ms))
                Nsteps_start = round(start/(self.dt/br2.ms))     
                
                fluctuant_array = np.zeros((Nsteps_total, self.network.basics.struct.Ncells))
                I_arr = np.asarray([0]*Nsteps_start+[function(t) for t in np.linspace(start, stop, Nsteps_total-Nsteps_start, endpoint=False)])
                
                if len(neuron_idcs)>0:
                    fluctuant_array[:, neuron_idcs] = repmat(I_arr, len(neuron_idcs), 1).transpose()
                        
            fluctuant_str = 'fluctuant_array(t, i)'
            
        else:
            fluctuant_array = None
            fluctuant_str = '0'
        
        self.fluctuant_array = fluctuant_array
        self.fluctuant_str = fluctuant_str
    
    def set_constant_stimuli(self):
        I_DC = np.zeros(self.network.basics.struct.Ncells_total)
        if self.constant_stimuli is not None:
            for target, value in self.constant_stimuli:
                I_DC[self.neuron_idcs(target)] = value
         
        self.neurons.I_DC = I_DC*br2.pA
    
    def get_membrane_events_dict(self):
        membrane_events = {}       
        for event in self.network.basics.equations.membr_events:        
            membrane_events[event] = self.network.basics.equations.membr_events[event]['condition']
        return membrane_events


    def set_membrane_events(self):
        self.event_monitors=NetworkHolder()
        
        for event in self.network.basics.equations.membr_events:
            self.event_monitors[event] = br2.EventMonitor(self.neurons, event, variables=self.network.basics.equations.membr_events[event]['vars'])
            self.neurons.run_on_event(event, self.network.basics.equations.membr_events[event]['reset'])    

    def set_params(self, obj, params, params_dict):
        for par in list(params_dict.keys()):
            
            unit = self.network.basics.membr.unit_dict[params_dict[par]['unit']]
            value = params_dict[par]['value']
            if isinstance(value, str):  
               setattr(obj, par, params.loc[dict(param=value)].values.astype(float)*unit)
            elif isinstance(value, int) or isinstance(value, float):
                 setattr(obj, par, value*unit)
                 
    def set_neuron_params(self):
        self.set_params(self.neurons, self.network.membr_params, self.network.basics.membr.name_units)
        
    def set_auxiliar_vars(self):
        self.neurons.I_ref = self.network.refractory_current.values*br2.pA     
        self.neurons.last_spike = -1000*br2.ms
        
    def set_initial_state(self):
        self.neurons.V = self.neurons.E_L
        self.neurons.w = 0 * br2.pA
    
    def set_STSP_vars(self):
        self.set_params(self.synapses, self.network.syn_params['STSP_params'], self.network.basics.syn.STSP.decl_vars)
       
    def set_syn_spike_params(self):
        self.set_params(self.synapses, self.network.syn_params['spiking'], self.network.basics.syn.spiking.names)
      
    def get_gsyn_amp(self):
        gsyn_amp={}
        for name in self.network.basics.syn.channels.names:
            factor = self.network.basics.syn.channels.gsyn_factor.loc[dict(channel=name, param=['factor'])].values
            tau_on, tau_off = self.network.basics.syn.channels.params.loc[dict(channel=name, param=['tau_on', 'tau_off'])].values            
            gsyn_amp[name] =  factor * tau_off * tau_on/(tau_off-tau_on)
        
        return gsyn_amp
    
    def set_channels(self):
        for name in self.network.basics.syn.channels.names:
            setattr(self.synapses, name, (self.network.syn_params['channel'].loc[dict(param='channel')].values==self.network.basics.syn.channels.names.index(name)).astype(int))
            self.synapses.gsyn_amp = self.synapses.gsyn_amp + getattr(self.synapses, name) * self.gsyn_amp[name]
            
    def set_delay(self):
        delay = self.network.syn_params['delay'].loc[dict(param='delay')].values.astype(float) * br2.ms
        
        for p in range(len(self.network.basics.equations.syn_pathway)):
            getattr(self.synapses, 'p{}'.format(p)).order = self.network.basics.equations.syn_pathway[p]['order']
            if self.network.basics.equations.syn_pathway[p]['delay']:
                getattr(self.synapses, 'p{}'.format(p)).delay = delay    
        
    def get_syn_pathway_dict(self):
        syn_dict_pathway = {}
        for p in range(len(self.network.basics.equations.syn_pathway)):
            syn_dict_pathway['p{}'.format(p)]=self.network.basics.equations.syn_pathway[p]['eq']
        return syn_dict_pathway 
    
    def get_ext_syn_pathway_dict(self):
        syn_dict_pathway = {}
        for p in range(len(self.network.basics.equations.ext_syn_pathway)):
            syn_dict_pathway['p{}'.format(p)]=self.network.basics.equations.ext_syn_pathway[p]['eq']
        return syn_dict_pathway 
    
    def set_synapses(self):
        
        self.gsyn_amp=self.get_gsyn_amp()   
        
        if self.network.syn_pairs.shape[1]>0:
            self.set_channels()         
            self.set_STSP_vars()
            self.set_syn_spike_params()
            self.set_delay()    
    
    def process_longrun(self, l, t1):
      
        if not os.path.isdir('longrun'):
            os.mkdir('longrun')
        name_units = self.network.basics.equations.var_units
        unit_dict = self.network.basics.membr.unit_dict
        
        longrun = self.longrun_monitor_control[l]
        variables = longrun['monitor'].keys()
        
        if ((longrun['stop'] is not None and longrun['start']<=t1<longrun['stop']) or (longrun['stop'] is None and longrun['start']<=t1)) and not longrun['longrun_monitor'].active:
      
            longrun['longrun_monitor'].active = True
            for var in variables:
                longrun['files'][var] = []
        elif longrun['longrun_monitor'].active:
          
            i = round((t1 - longrun['start'])/longrun['interval'])
            for var in variables:
                mon_var = getattr(longrun['longrun_monitor'], var)/unit_dict[name_units[var]]
                if longrun['population_agroupate'] is not None and var!='t':
                    if longrun['population_agroupate'] == 'mean':
                        mon_var = np.mean(mon_var, axis=0)
                    elif longrun['population_agroupate'] == 'sum':
                        mon_var = np.sum(mon_var, axis=0)
          
            
                with open('longrun/{}_{}_{}.npy'.format(var, i, l),'wb') as f:
                    np.save(f, mon_var)
                longrun['files'][var].append('longrun/{}_{}_{}.npy'.format(var, i, l))
            idc = longrun['longrun_monitor'].record
            v = longrun['longrun_monitor'].record_variables
            source = longrun['longrun_monitor'].source
            
            self.net.remove(longrun['longrun_monitor'])
            longrun['longrun_monitor'] = br2.StateMonitor(source, v, idc, dt=self.dt)
            self.net.add(longrun['longrun_monitor'])
        
        
        if longrun['stop'] is not None and t1 >= longrun['stop'] and  longrun['longrun_monitor'].active:
          
            longrun['longrun_monitor'].active=False
            self.net.remove(longrun['longrun_monitor'])
    
    def restore_longrun(self, erase=True):
        
        name_units = self.network.basics.equations.var_units
        unit_dict = self.network.basics.membr.unit_dict
        for longrun in self.longrun_monitor_control:
            for var in list(longrun['files'].keys()):
                longvar = []
                for file in longrun['files'][var]:                                   
                    longvar.append(np.load(file))  
                
                longrun['monitor'][var] = np.concatenate(longvar) if longrun['population_agroupate'] is not None  or var=='t' else  np.concatenate(longvar, axis=1)
                longrun['monitor'][var] = longrun['monitor'][var]*unit_dict[name_units[var]]
                
        if len(self.longrun_monitor_control) and erase:
            shutil.rmtree('longrun')
    
    def set_custom_stimuli(self, name, Nsource, channel, spike_idcs, spike_times, pairs_connected, gmax, pfail):
        
        if isinstance(channel, str):
            channel = [channel]        
        
        Nconnected = pairs_connected.shape[1]
        
        channel_arr = []
        for ch in channel:
            channel_arr.append([ch]*Nconnected)
        channel_arr = np.concatenate(channel_arr)
        
        pairs_connected = repmat(pairs_connected, 1, len(channel))
        
        # gmax_arr = []
        if isinstance(gmax, (int, float)):
            gmax = [gmax]
        
        if len(gmax)==1:
            gmax_arr = gmax*Nconnected*len(channel)
        elif len(gmax)==Nconnected:   
            gmax_arr = gmax*len(channel)
            
        if isinstance(pfail, (int, float)):
            pfail = [pfail]
            
        if len(pfail)==1:
            pfail_arr = pfail*Nconnected*len(channel)
        elif len(pfail)==Nconnected:   
            pfail_arr = pfail*len(channel)
        
            
        
        self.external_stimuli[name] = NetworkHolder()
        self.external_stimuli[name].generator = br2.SpikeGeneratorGroup(Nsource, spike_idcs, spike_times*br2.ms, dt=self.dt)
        ext_syn_model = self.network.basics.equations.ext_syn_model
        ext_syn_dict = self.get_ext_syn_pathway_dict()
        
        self.external_stimuli[name].synapses = br2.Synapses(self.external_stimuli[name].generator, self.neurons, model=ext_syn_model, on_pre=ext_syn_dict, dt=self.dt)
        self.external_stimuli[name].synapses.connect(i=pairs_connected[1,:], j=pairs_connected[0,:])

        for ch_name in self.network.basics.syn.channels.names:
            setattr(self.external_stimuli[name].synapses, ch_name, (channel_arr==ch_name).astype(int))
            self.external_stimuli[name].synapses.gsyn_amp = self.external_stimuli[name].synapses.gsyn_amp + getattr(self.external_stimuli[name].synapses, ch_name) * self.gsyn_amp[ch_name]
        
        self.external_stimuli[name].synapses.pfail = pfail_arr
        self.external_stimuli[name].synapses.gmax = gmax_arr *br2.nS 
        self.external_stimuli[name].spikemonitor = br2.SpikeMonitor(self.external_stimuli[name].generator)
        
        self.net.add(self.external_stimuli[name].generator, self.external_stimuli[name].synapses, self.external_stimuli[name].spikemonitor)
 
    def set_poisson_stimuli(self, stimulator_name, Nsource, channel, target_idc, pCon, rate, start, stop, gmax, pfail):
    
        spike_times = []
        spike_idcs = []
     
        source_idc = np.arange(Nsource)
        Ntarget = len(target_idc)
        for src in range(Nsource):
            lambd = 1000/rate
            upper_est = int(round(2*(stop-start)/lambd, 0))
            spike_array =  start + np.cumsum(np.clip(-lambd *np.log(1-np.random.rand(upper_est)), self.dt/br2.ms, np.NaN))
            
            spikeact = spike_array[spike_array < stop]
                    
            spike_times.extend(spikeact)
            spike_idcs.extend([src]*len(spikeact))
        spike_idcs = np.asarray(spike_idcs)
        
        Nconnections = int(round(pCon*Ntarget*Nsource, 0))
        
        local_pairs = np.arange(Nsource*Ntarget)
        np.random.shuffle(local_pairs)
        local_pairs = local_pairs[:Nconnections]
        local_target_idc = local_pairs//Nsource
        local_source_idc = local_pairs % Nsource
        
        pairs_connected = np.zeros((2,Nconnections))
        pairs_connected[0,:] = target_idc[local_target_idc]
        pairs_connected[1,:] = source_idc[local_source_idc]
        pairs_connected=pairs_connected.astype(int)
        
        
        # gmax_arr = np.array([gmax]*Nconnections)
        # pfail_arr = np.array([pfail]*Nconnections)
        
        # if isinstance(channel, str):
        #     channel = [channel]
         
        # pairs_connected = repmat(pairs_connected, 1, len(channel))
        # gmax_arr = repmat(gmax_arr, 1, len(channel))
        # pfail_arr = repmat(pfail_arr, 1, len(channel))
        
        # channel_arr = []
        # for ch in channel:
        #     channel_arr.append([ch]*Nconnections)
        # channel_arr = np.concatenate(channel_arr)
        
        
        # stimulator_setup = StimulatorSetup(Nsource, spike_idcs, spike_times, pairs_connected, channel, gmax, pfail)
        
        self.set_custom_stimuli(stimulator_name, Nsource, channel, spike_idcs, spike_times, pairs_connected, gmax, pfail)
    
    
    def set_regular_stimuli(self, stimulator_name, Nsource, channel, target_idc, pCon, rate, start, stop, gmax, pfail):

        source_idc = np.arange(Nsource)
        Nspikes=round((stop-start)*rate/1000)

        if 1000/rate < self.dt/br2.ms:
            raise ValueError('ERROR: time step is larger than spikes intervals')

        spike_times = list(np.linspace(start, stop, Nspikes, endpoint=False))*Nsource
        spike_idcs=np.concatenate(([[src]*Nspikes for src in source_idc])).astype(int)

        Nconnections = int(round(pCon*len(target_idc)*Nsource, 0))

        local_pairs = np.arange(Nsource*len(target_idc))
        np.random.shuffle(local_pairs)
        local_pairs = local_pairs[:Nconnections]
        local_target_idc = local_pairs//Nsource
        local_source_idc = local_pairs % Nsource

        pairs_connected = np.zeros((2,Nconnections))
        pairs_connected[0,:] = target_idc[local_target_idc]
        pairs_connected[1,:] = source_idc[local_source_idc]
        pairs_connected=pairs_connected.astype(int)

        # stimulator_setup = StimulatorSetup(Nsource, spike_idcs, spike_times, pairs_connected, channel, gmax, pfail)
        # return stimulator_setup
        self.set_custom_stimuli(stimulator_name, Nsource, channel, spike_idcs, spike_times, pairs_connected, gmax, pfail)
    
    def set_monitors(self, monitor, variables, start, stop):
        
        if self.transient>0:
            monitor.active = False
        for var in variables:
            self.recorded[var] = monitor
            
        if start is None:
            start = self.transient
        
        if start in self.monitor_schedule:
            self.monitor_schedule[start]['start'].append(monitor)
        else:
            self.monitor_schedule[start] = dict(start=[monitor], stop=[])
            
        if stop is not None:
            if stop in self.monitor_schedule:
                self.monitor_schedule[stop]['stop'].append(monitor)
            else:
                self.monitor_schedule[stop] = dict(start=[], stop=[monitor])
                
        self.net.add(monitor)
 
    def set_neuron_monitors(self, name, variables, groupstripe_list, start=None, stop=None):
        
        if isinstance(variables, str):
            variables = [variables]
        
        neuron_idc =self.neuron_idcs(groupstripe_list)
              
        self.neuron_monitors[name] =  br2.StateMonitor(self.neurons, variables, neuron_idc, dt=self.dt)
        self.set_monitors(self.neuron_monitors[name], variables, start, stop)
        
    def set_synapse_monitors(self, name, variables, target_groupstripe_list, source_groupstripe_list, start=None, stop=None):
      
        if isinstance(variables, str):
            variables = [variables]
        
        syn_idc = self.syn_idcs_from_groups(target_groupstripe_list, source_groupstripe_list)
          
        self.synapse_monitors[name] =  br2.StateMonitor(self.neurons, variables, syn_idc, dt=self.dt)
        self.set_monitors(self.synapse_monitors[name], variables, start, stop)      
    
    def set_longrun_monitors(self, monitor, longrun_monitor, variables, interval, start, stop, population_agroupate):
        
        if self.transient>0:
            longrun_monitor.active = False
        for var in variables:
            self.recorded[var] = monitor
            
        if start is None:
            start = self.transient
        
        self.longrun_monitor_control.append(dict(start=start, stop=stop, interval=interval,
                                       longrun_monitor=longrun_monitor, monitor=monitor,
                                       population_agroupate=population_agroupate, files={}))
        
        
        self.net.add(longrun_monitor)    
    
    def set_longrun_neuron_monitors(self, name, variables, groupstripe_list,  interval, start=None, stop=None, population_agroupate=None):
        
        if isinstance(variables, str):
            variables = [variables]
        
        neuron_idc =self.neuron_idcs(groupstripe_list)
              
        self.longrun_neuron_monitors[name] =  br2.StateMonitor(self.neurons, variables, neuron_idc, dt=self.dt)
        self.longrun_neuron_monitors[name].active = False
        self.neuron_monitors[name] = NetworkHolder()
        self.neuron_monitors[name].t = None
        for var in variables:
            self.neuron_monitors[name][var] = None
            
        self.set_longrun_monitors(self.neuron_monitors[name], self.longrun_neuron_monitors[name],
                                  variables, interval, start, stop, population_agroupate)
        
    def set_longrun_synapse_monitors(self, name, variables, target_groupstripe_list, source_groupstripe_list,  interval, start=None, stop=None, population_agroupate=None):
        
        if isinstance(variables, str):
            variables = [variables]
        
        syn_idc = self.syn_idcs_from_groups(target_groupstripe_list, source_groupstripe_list)
          
        self.longrun_synapse_monitors[name] =  br2.StateMonitor(self.neurons, variables, syn_idc, dt=self.dt)
        self.longrun_synapse_monitors[name].active = False
        self.synapse_monitors[name] = NetworkHolder()
        self.synapse_monitors[name].t = None
        for var in variables:
            self.synapse_monitors[name][var] = None
            
        self.set_longrun_monitors(self.synapse_monitors[name], self.longrun_synapse_monitors[name],
                                  variables, interval, start, stop, population_agroupate)
        
    def group_idcs(self, group):
        return self.network.group_idcs(group)

    def neuron_idcs(self, groupstripe_list):    
        return self.network.neuron_idcs(groupstripe_list)
    
    def syn_idcs_from_neurons(self, target, source, channel=None):
        return self.network.syn_idcs_from_neurons(target, source, channel)
    
    def syn_idcs_from_groups(self, target_groupstripe_list, source_groupstripe_list, channel=None):
        return self.network.syn_idcs_from_groups(target_groupstripe_list, source_groupstripe_list, channel)
        
    

    
    def spiking_idcs(self, comparisons, groupstripe=None, tlim=None):
        
        if isinstance(comparisons[0], np.ufunc):
            comparisons = [comparisons]
            
        rate = self.spiking_rate(tlim)
        
        idc_bool_list = []
        for relation, comparison_rate in comparisons: 
            idc_bool_list.append(relation(rate, comparison_rate))
        
        sp_idc = np.arange(self.neurons.N)
        
        if groupstripe is not None:
            neuronidc = self.neuron_idcs(groupstripe)
            sp_isin = np.isin(sp_idc, neuronidc)
            idc_bool_list.append(sp_isin)
            
        idc_bool = np.isfinite(sp_idc)
        for item in idc_bool_list:
            idc_bool = idc_bool & item
            
        sp_idc= sp_idc[idc_bool]
        
        return sp_idc
   
    def spiking_count(self, neuron_idcs=None, tlim=None, delta_t=None):
        
        if tlim is None and delta_t is None:
            count = self.spikemonitor.count
            if neuron_idcs is not None:
                count = count[neuron_idcs]
            return count
        else:
            if tlim is None:
                t0, t1 = cortex.transient, cortex.neurons.t/br2.ms
            else:
                t0,t1 = tlim
            
            if neuron_idcs is None:
                neuron_idcs = np.arange(self.neurons.N)
                
            spike_trains = self.spikemonitor.spike_trains()
       
            if delta_t is None:
                count=[]
                for idc in neuron_idcs:
                    spiking = spike_trains[idc]/br2.ms
                    count.append(len(spiking[(spiking>=t0)&(spiking<t1)]))
            else:
                N_intervals = np.ceil((t1-t0)/delta_t).astype(int)
                count = [[] for i in  range(len(neuron_idcs))]
                for i in range(len(neuron_idcs)):
                    spiking = spike_trains[neuron_idcs[i]]/br2.ms                
                    for n in range(N_intervals):
                        count[i].append(len(spiking[(spiking>=t0+n*delta_t)&(spiking<t0+(n+1)*delta_t)]))
                
            return np.asarray(count)
    
    
    def spiking_rate(self, neuron_idcs=None, tlim=None):
        count = self.spiking_count(neuron_idcs, tlim)
        if tlim is None:
            return count/(self.neurons.t/br2.second - self.transient/1000)
        else:
            t0,t1 = tlim
            return count/((t1-t0)/1000)
       
    def rheobase(self, neuron_idcs=None):     
        return self.network.rheobase(neuron_idcs)

     
@dataclass
class NetworkHolder(BaseClass):
    pass
    
    
@dataclass
class VarHolder(NetworkHolder):
    pass
    
    def var(self, var):
        return getattr(self[var], var)
    
    def t(self, var):
        return self[var].t
    
        
@dataclass
class StimulatorSetup(BaseClass):
    Nsource: int
    spike_idcs: np.array
    spike_times: np.array
    pairs: np.array
    channel:np.array
    gmax: np.array
    pfail: np.array
    
    def __post_init__(self):
        self.channel = np.asarray(self.channel)
        self.spike_idcs = np.asarray(self.spike_idcs)
        self.spike_times = np.asarray(self.spike_times)
        self.pfail= np.asarray(self.pfail)
        
    def __add__(self, instance2):
        Nsource = self.Nsource + instance2.Nsource
        spike_idcs = np.concatenate((self.spike_idcs, instance2.spike_idcs))
        spike_times = np.concatenate((self.spike_times, instance2.spike_times))
        pairs = np.concatenate((self.pairs, instance2.pairs))
        channel = np.concatenate((self.channel, instance2.channel))
        gmax=np.concatenate((self.gmax, instance2.gmax))
        pfail=np.concatenate((self.pfail, instance2.pfail))
        return StimulatorSetup(Nsource, spike_idcs, spike_times, pairs, channel, gmax, pfail)

    



    
@time_report('Task')
def run():
    
    duration=2000
    seed = None

    transient=0
  
    constant_stimuli = [
                        [('PC', 0), 250],
                        [('IN', 0), 200]
                        ]
    
    
    basics_scales = {'membr_param_std': [(dict(param=MEMB_PARAM_NAMES, group=SETS['ALL']), 0)] }

    # fluctuant_stimuli = [[('PC_L23', 0), lambda t: 50+50*np.sin(1*2*np.pi*(t/1000)) , 600, 1200]]
                          
    # cortex = Cortex(network, constant_stimuli, method=method, dt=dt, seed=seed)
    # cortex.set_neuron_monitors('I_tot', 'I_tot', ('ALL',0))
    cortex=Cortex.setup(Ncells=1000, Nstripes=1, constant_stimuli=constant_stimuli, method='rk4',
                        dt=0.05, transient=1000,seed=seed, basics_scales=basics_scales)
    # cortex = Cortex.load(path='Network2', constant_stimuli=constant_stimuli, method='rk4', dt=0.05, transient=transient)
    cortex.run(duration)
    # I_inj_mon=cortex.neuron_monitors['I_tot']
    # plt.plot(I_inj_mon.t/br2.ms, np.mean(I_inj_mon.I_tot, axis=0)/br2.pA)
    return cortex
    # tools.raster_plot(cortex)
    
    # elapsed = tm.time() - start_time
    # print("--- {} min {:.2f} seconds ---".format(elapsed//60, elapsed%60))
    
if __name__ == '__main__':
    
    cortex = run()    