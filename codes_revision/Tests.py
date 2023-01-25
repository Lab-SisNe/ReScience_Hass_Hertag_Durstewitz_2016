from CortexSetup import Cortex
from Auxiliar import set_simulation_dir, time_report
from NetworkSetup import network_setup
from BasicsSetup import GROUPS_SETS as SETS, MEMB_PARAM_NAMES
import AnalysesTools as Tools
import numpy as np
import Figures
import os

# 10 min and 35.42 s.
@time_report()
def task1(simulation_dir, seed=None):
    
    Ncells=1000
    Nstripes=1
    constant_stimuli = [
                        [('PC', 0), 250],
                        [('IN', 0), 200]
                        ]
    method='rk4'
    dt=0.05
    transient=1000   
    seed = seed
    
    duration=31000
                          
    cortex = Cortex.setup(Ncells=Ncells, Nstripes=Nstripes, constant_stimuli=constant_stimuli,method=method, dt=dt,
                          transient=transient, seed=seed)
    cortex.set_longrun_neuron_monitors('I_tot', 'I_tot', ('ALL',0),  5000, start=1000, stop=31000, population_agroupate='sum')
    cortex.set_longrun_neuron_monitors('V', 'V', ('ALL',0),  5000, start=1000, stop=31000)
    
    cortex.run(duration)
    
    Figures.fig01(cortex, (duration-6000, duration), simulation_dir)
    
    aliases = ['PC', 'fast-spiking cells','bitufted cells', 'basket cells', 'Matinotti cells']
    Tools.get_membr_params(cortex,[('PC',0), ('IN_L_both',0), ('IN_CL_both',0), ('IN_CC',0), ('IN_F',0)], alias_list = aliases, file=simulation_dir+'//membr_params.txt');
    Tools.get_spiking(cortex, 0.33, ('PC', 0), file=simulation_dir+'\\spiking.txt');
    
    Tools.comp_membrparam_rategroup(cortex, 0.33, [('PC_L23',0), ('PC_L5', 0), ('PC', 0), ('ALL',0)], file=simulation_dir+'\\membr_params_comparison.txt');
    for channel in cortex.network.basics.syn.channels.names:
        Tools.contingency(cortex, 0.33, [('PC_L23',0), ('PC_L5', 0), ('PC', 0), ('ALL',0)], ('ALL', 0), channel=channel,file=simulation_dir+'\\pcon_contingency_{}.txt'.format(channel) );
        Tools.comp_synparam_rategroup(cortex, 0.33, [('PC_L23',0), ('PC_L5', 0), ('PC', 0), ('ALL',0)], ('ALL', 0), channel=channel, file=simulation_dir+'\\syn_params_comparison_{}.txt'.format(channel));
        
    
    spikingPC = cortex.spiking_idcs((np.greater_equal, 0.33), ('PC',0))
    spikingALL= cortex.spiking_idcs((np.greater_equal, 0.33), ('ALL', 0))
    
    lags = np.arange(50)
    tlim=(transient, transient+15000)
    
    C_lags, autoC_mean, autoC_std, crossC_mean, crossC_std= Tools.get_correlations(cortex, idcs=spikingALL, tlim=tlim, delta_t=2, lags=lags, file=simulation_dir+'\\Correlations.txt', display=True, display_interval=5)
    ISImean, ISICV = Tools.get_ISI_stats(cortex, neuron_idcs=spikingALL, tlim=tlim,  savetxt=simulation_dir+'\\ISI_stats.txt')
    Correlation_sigma = 1
    Figures.fig02(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, simulation_dir)
    Figures.fig03(C_lags, autoC_mean, simulation_dir)
    
    fq, pwr = Tools.get_LFP_SPD(cortex, log=True, sigma=2)
    Figures.fig04(fq, pwr, simulation_dir)
    
    if not os.path.isdir(simulation_dir+'\\Reports'):
        os.mkdir(simulation_dir+'\\Reports')
    
    VstdALL = Tools.get_V_stats(cortex, spikingALL, file=simulation_dir + '\\Reports\\VALLstats.txt')[1]
    VstdPC =Tools.get_V_stats(cortex, spikingPC, file=simulation_dir+'\\Reports\\VPCstats.txt')[1]
    Figures.fig05(VstdALL, VstdPC, simulation_dir)
# 81 min and 12.52 s.
@time_report()
def task2(simulation_dir):
    
    with open('Original_spiketime.txt', 'r') as f:
        text = f.read()
    
    
    spike_trains = [np.asarray(train.split(',')).astype(float) for train in text.split(';')]
    
    sp_train =[]
    tlim=(1000,16000)
    t0, t1 = tlim
    minfq = 0.33
    mincount = (t1-t0)/1000 * minfq
    for train in spike_trains:   
        if len(train[(train>=t0) & (train<t1)])>=mincount:
            sp_train.append(train)
    spike_trains = sp_train   

    
    tlim=(1000,16000)
    
    if not os.path.isdir(simulation_dir+'\\Reports'):
        os.mkdir(simulation_dir+'\\Reports')
    ISImean, ISICV = Tools.get_ISI_stats_from_spike_trains(spike_trains, tlim=tlim,  savetxt=simulation_dir+'\\Reports\\original_ISI.txt')
    
    lags= np.arange(50)
    
    
    C_lags, autoC_mean, autoC_std,crossC_mean, crossC_std = Tools.get_correlations_from_spike_trains(spike_trains, tlim=tlim, delta_t=2, lags=lags,file=simulation_dir+'\\Reports\\Original_spikes.txt', display=True, display_interval=5)  
    
    
    Correlation_sigma = 1
    Figures.fig06(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, simulation_dir)
    Figures.fig07(C_lags, autoC_mean, simulation_dir)
    
    
    time = np.load('Original_t.npy')
    Itot = np.load('Original_I.npy')
    
    time_bins = np.floor(time*20).astype(int)
    time_bins, unique_idc = np.unique(time_bins, return_index=True)
    
    Itot = Itot[unique_idc]
    frequency = 1000/0.05
    fq, pwr = Tools.get_LFP_SPD_from_Itotarray(Itot, frequency, log=True, sigma=0)
    
    Figures.fig08(fq, pwr, simulation_dir)

@time_report()
def task3(simulation_dir, seed=None):

    Ncells=1000
    Nstripes=1
    constant_stimuli = [
                        [('PC', 0), 250],
                        [('IN', 0), 200]
                        ]
    method='gsl_rk2'
    dt=0.05
    transient=1000   
    seed = seed
    
    duration=16000
                          
    cortex = Cortex.setup(Ncells=Ncells, Nstripes=Nstripes, constant_stimuli=constant_stimuli,method=method, dt=dt,
                          transient=transient, seed=seed)
    # cortex.save('networkAAA')
    cortex.set_longrun_neuron_monitors('I_tot', 'I_tot', ('ALL',0),  5000, start=1000, stop=31000, population_agroupate='sum')
    cortex.set_longrun_neuron_monitors('V', 'V', ('ALL',0),  5000, start=1000, stop=31000)

    cortex.run(duration)
    
    spikingALL= cortex.spiking_idcs((np.greater_equal, 0.33), ('ALL', 0))

    Figures.fig09(cortex, (duration-6000, duration), simulation_dir)
    
    lags = np.arange(50)
    tlim=(transient, transient+15000)
    if not os.path.isdir(simulation_dir+'\\Reports'):
        os.mkdir(simulation_dir+'\\Reports')
    C_lags, autoC_mean, autoC_std, crossC_mean, crossC_std= Tools.get_correlations(cortex, idcs=spikingALL, tlim=tlim, delta_t=2, lags=lags, file=simulation_dir+'\\Reports\\Correlationsrk2.txt', display=True, display_interval=5)
    ISImean, ISICV = Tools.get_ISI_stats(cortex, neuron_idcs=spikingALL, tlim=tlim, savetxt=simulation_dir+'\\Reports\\ISIrk2_stats.txt')
    Correlation_sigma = 1
    
    Figures.fig10(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, simulation_dir)
    Figures.fig11(C_lags, autoC_mean, simulation_dir)
    
    fq, pwr = Tools.get_LFP_SPD(cortex, log=True, sigma=2)
    Figures.fig12(fq, pwr, simulation_dir)
    
    
@time_report()
def task4(simulation_dir, seed=None):
    
    Ncells=1000
    Nstripes=1
    constant_stimuli = [
                        [('PC', 0), 250],
                        [('IN', 0), 200]
                        ]
    method='rk4'
    dt=0.01
    transient=1000   
    seed = seed
    
    duration=2500
                          
    cortex = Cortex.setup(Ncells=Ncells, Nstripes=Nstripes, constant_stimuli=constant_stimuli,method=method, dt=dt,
                          transient=transient, seed=seed)
     
    cortex.set_neuron_monitors('I_tot', 'I_tot', ('ALL',0))
    PC_L23 = cortex.neuron_idcs(('PC_L23',0))
    
    pCon_reg = 0.1
    pulse0, rate0 = (1100, 1105), 50000
    pulse1, rate1 = (2100, 2105), 100000
    gmax_reg = 0.1
    pfail_reg=0
    cortex.set_regular_stimuli('regular1', 1, ['AMPA', 'NMDA'], PC_L23, pCon=pCon_reg, rate=rate0, start=pulse0[0], stop=pulse0[1], gmax=gmax_reg, pfail=pfail_reg)
    cortex.set_regular_stimuli('regular2', 1, ['AMPA', 'NMDA'], PC_L23, pCon=pCon_reg, rate=rate1, start=pulse1[0], stop=pulse1[1], gmax=gmax_reg, pfail=pfail_reg)
    cortex.run(duration)
    
    Figures.fig13(cortex, pulse0, pulse1, simulation_dir)
    
    PC_L23 = cortex.neuron_idcs(('PC_L23', 0))
    PC_L5 = cortex.neuron_idcs(('PC_L5', 0))
    
    delta_t=1
    
    tpulse0 = np.arange(pulse0[0], pulse0[1]+50, delta_t)
    PCL23spikes_pulse0 = np.sum(cortex.spiking_count(neuron_idcs=PC_L23,tlim=(pulse0[0], pulse0[0]+50), delta_t=delta_t), axis=0)
    PCL5spikes_pulse0 = np.sum(cortex.spiking_count(neuron_idcs=PC_L5,tlim=(pulse0[0], pulse0[0]+50), delta_t=delta_t), axis=0)
    
    tpulse1 = np.arange(pulse1[0], pulse1[1]+50, delta_t)
    PCL23spikes_pulse1 = np.sum(cortex.spiking_count(neuron_idcs=PC_L23,tlim=(pulse1[0], pulse1[0]+50), delta_t=delta_t), axis=0)
    PCL5spikes_pulse1 = np.sum(cortex.spiking_count(neuron_idcs=PC_L5,tlim=(pulse1[0], pulse1[0]+50), delta_t=delta_t), axis=0)
    if not os.path.isdir(simulation_dir+'\\Reports'):
        os.mkdir(simulation_dir+'\\Reports')
        
    with open(simulation_dir+'\\Reports\\Regularpulses_spikingcounts.txt', 'w') as f:
        print('Regular pulses on PC_L23', end='\n\n', file=f)
        print('gmax:', gmax_reg, file=f)
        print('pCon:', pCon_reg, file=f)
        print('pfail', pfail_reg, file=f, end='\n\n')
        
        print('Pulse0', file=f)
        print('Start: {} ms, stop: {} ms'.format(*pulse0), file=f)
        print('Rate: {} Hz'.format(rate0),file=f, end='\n\n')
        
        print('Spike count on PC_L23', file=f)
        time_spikes = list(zip(tpulse0, PCL23spikes_pulse0))
        for i in range(len(time_spikes)):
            print('t: {} ms, spikes: {}'.format(*time_spikes[i]), file=f)
        print('Total spikes on PC_L23: {}'.format(np.sum(PCL23spikes_pulse0)), end='\n\n', file=f)
        
        print('Spike count on PC_L5', file=f)
        time_spikes = list(zip(tpulse0, PCL5spikes_pulse0))
        for i in range(len(time_spikes)):
            print('t: {} ms, spikes: {}'.format(*time_spikes[i]), file=f)
        print('Total spikes on PC_L5: {}'.format(np.sum(PCL5spikes_pulse0)), end='\n\n', file=f)
        
        
        print('Pulse1', file=f)
        print('Start: {} ms, stop: {} ms'.format(*pulse1), file=f)
        print('Rate: {} Hz'.format(rate1),file=f, end='\n\n')
        
        print('Spike count on PC_L23', file=f)
        time_spikes = list(zip(tpulse1, PCL23spikes_pulse1))
        for i in range(len(time_spikes)):
            print('t: {} ms, spikes: {}'.format(*time_spikes[i]), file=f)
        print('Total spikes on PC_L23: {}'.format(np.sum(PCL23spikes_pulse1)), end='\n\n', file=f)
        
        print('Spike count on PC_L5', file=f)
        time_spikes = list(zip(tpulse1, PCL5spikes_pulse1))
        for i in range(len(time_spikes)):
            print('t: {} ms, spikes: {}'.format(*time_spikes[i]), file=f)
        print('Total spikes on PC_L5: {}'.format(np.sum(PCL5spikes_pulse1)), end='\n\n', file=f)
        
@time_report()
def task5(simulation_dir, seed=None):
    basics_scales = {'membr_param_std': [(dict(param=MEMB_PARAM_NAMES, group=SETS['ALL']), 0.8)] }

    Ncells=1000
    Nstripes=1
    constant_stimuli = [
                        [('PC', 0), 250],
                        [('IN', 0), 200]
                        ]
    method='rk4'
    dt=0.01
    transient=1000   
    seed = seed
                          
    cortex = Cortex.setup(Ncells=Ncells, Nstripes=Nstripes, constant_stimuli=constant_stimuli,method=method, dt=dt,
                          basics_scales=basics_scales, transient=transient, seed=seed)
    
    
    duration=1200 
    pCon_reg = 0.1
    pulse=(1100, 1105)
    rate=100000
    gmax_reg = 0.1
    pfail_reg=0
    
    PC_L23 = cortex.neuron_idcs(('PC_L23',0))
    PC_L5 = cortex.neuron_idcs(('PC_L5',0))
    
    cortex.set_regular_stimuli('regular1', 1, ['AMPA', 'NMDA'], PC_L23, pCon=pCon_reg, rate=rate, start=pulse[0], stop=pulse[1], gmax=gmax_reg, pfail=pfail_reg)
    
    cortex.run(duration)
    
    Figures.fig14(cortex, pulse, simulation_dir)
    
    delta_t=1
    
    tpulse = np.arange(pulse[0], pulse[0]+50, delta_t)
    PCL23spikes_pulse = np.sum(cortex.spiking_count(neuron_idcs=PC_L23,tlim=(pulse[0], pulse[0]+50), delta_t=delta_t), axis=0)
    PCL5spikes_pulse = np.sum(cortex.spiking_count(neuron_idcs=PC_L5,tlim=(pulse[0], pulse[0]+50), delta_t=delta_t), axis=0)
    
    if not os.path.isdir(simulation_dir+'\\Reports'):
        os.mkdir(simulation_dir+'\\Reports')
        
    with open(simulation_dir+'\\Reports\\Regularpulses_reducedstd__spikingcounts.txt', 'w') as f:
        print('Regular pulses on PC_L23', end='\n\n', file=f)
        print('gmax:', gmax_reg, file=f)
        print('pCon:', pCon_reg, file=f)
        print('pfail', pfail_reg, file=f, end='\n\n')
        
        print('Pulse', file=f)
        print('Start: {} ms, stop: {} ms'.format(*pulse), file=f)
        print('Rate: {} Hz'.format(rate),file=f, end='\n\n')
        
        print('Spike count on PC_L23', file=f)
        time_spikes = list(zip(tpulse, PCL23spikes_pulse))
        for i in range(len(time_spikes)):
            print('t: {} ms, spikes: {}'.format(*time_spikes[i]), file=f)
        print('Total spikes on PC_L23: {}'.format(np.sum(PCL23spikes_pulse)), end='\n\n', file=f)
        
        print('Spike count on PC_L5', file=f)
        time_spikes = list(zip(tpulse, PCL5spikes_pulse))
        for i in range(len(time_spikes)):
            print('t: {} ms, spikes: {}'.format(*time_spikes[i]), file=f)
        print('Total spikes on PC_L5: {}'.format(np.sum(PCL5spikes_pulse)), end='\n\n', file=f)

@time_report()
def task6(simulation_dir, seed=None):

    Ntrials = 5
    Membr_std_list = np.asarray([0,20,40,60,80,100])/100
    
    
    PC_L23_results = [[] for i in range(len(Membr_std_list))]
    PC_L5_results = [[] for i in range(len(Membr_std_list))]
    
    
    
    constant_stimuli = [
                        [('PC', 0), 250],
                        [('IN', 0), 200]
                        ]
    
    method = 'rk4'
    dt = 0.01
    transient= 1000
    duration = 1200
    
    

    Ncells=1000
    Nstripes=1
    constant_stimuli = [
                        [('PC', 0), 250],
                        [('IN', 0), 200]
                        ]

    
    for trial in range(Ntrials):
        for std in range(len(Membr_std_list)):
            print_out1a = '||   REPORT: Trial ' + str(trial)
            print_out2a = '||   REPORT: std_scale ' + str(Membr_std_list[std])
            print_out1 = print_out1a + ' '*(max(len(print_out1a), len(print_out2a)) - len(print_out1a)) + '   ||'
            print_out2 = print_out2a + ' '*(max(len(print_out1a), len(print_out2a)) - len(print_out2a)) + '   ||'
            print("="*len(print_out1))
            print(print_out1)
            print(print_out2)
            print("="*len(print_out1))
            
            basics_scales = {'membr_param_std': [(dict(param=MEMB_PARAM_NAMES, group=SETS['ALL']), Membr_std_list[std])] }



            cortex = Cortex.setup(Ncells=Ncells, Nstripes=Nstripes, constant_stimuli=constant_stimuli,method=method, dt=dt,
                                  basics_scales=basics_scales, transient=transient, seed=seed)
            
            
            
            PC_L23 = cortex.neuron_idcs(('PC_L23',0))
            PC_L5 = cortex.neuron_idcs(('PC_L5',0))
            
            pCon_reg = 0.1
            pulse=(1100, 1105)
            rate=100000
            gmax_reg = 0.1
            pfail_reg=0
            
            cortex.set_regular_stimuli('regular', 1, ['AMPA', 'NMDA'], PC_L23, pCon=pCon_reg, rate=rate, start=pulse[0], stop=pulse[1], gmax=gmax_reg, pfail=pfail_reg)
 
            cortex.run(duration)        
 
            
            
            
            
            PC_L23_results[std].append(np.sum(cortex.spiking_count(neuron_idcs=PC_L23,tlim=(pulse[0], pulse[0]+50))))
            PC_L5_results[std].append(np.sum(cortex.spiking_count(neuron_idcs=PC_L5,tlim=(pulse[0], pulse[0]+50))))

    Figures.fig15(Membr_std_list, PC_L23_results, PC_L5_results, simulation_dir)

@time_report()
def task7(simulation_dir, seed=None):

    Ncells=1000
    Nstripes=1
    constant_stimuli = [
                        [('PC', 0), 250],
                        [('IN', 0), 200]
                        ]
    method='rk4'
    dt=0.05
    transient=1000   
    seed = seed
    
    duration=3000
    
    cortex = Cortex.setup(Ncells=Ncells, Nstripes=Nstripes, constant_stimuli=constant_stimuli,method=method, dt=dt,
                          transient=transient, seed=seed)
     
    
    pCon_reg = 0.1
    pulse0=(1100, 1200)
    pulse1=(2100, 2200)
    rate=30
    gmax_reg = 0.1
    pfail_reg=0
    PC_L23 = cortex.neuron_idcs(('PC_L23',0))
    PC_L5 = cortex.neuron_idcs(('PC_L5',0))
    
    
    cortex.set_regular_stimuli('poisson1', 100, ['AMPA', 'NMDA'], PC_L23, pCon=pCon_reg, rate=rate, start=pulse0[0], stop=pulse0[1], gmax=gmax_reg, pfail=pfail_reg)
    cortex.set_regular_stimuli('poisson1', 100, ['AMPA', 'NMDA'], PC_L5, pCon=pCon_reg, rate=rate, start=pulse1[0], stop=pulse1[1], gmax=gmax_reg, pfail=pfail_reg)
    
    cortex.run(duration)
    
    Figures.fig16(cortex, pulse0, pulse1, simulation_dir)

@time_report()
def task8(simulation_dir, seed=None):
    
    basics_scales = {'gmax_mean': [(dict(target=SETS['ALL'], source=SETS['IN']), 0.4)]  }
    
    
    Ncells=1000
    Nstripes=1
    constant_stimuli = [
                        [('PC', 0), 250],
                        [('IN', 0), 200]
                        ]
    method='rk4'
    dt=0.05
    transient=1000   
    seed = seed
    
    duration=3000
    
    cortex = Cortex.setup(Ncells=Ncells, Nstripes=Nstripes, constant_stimuli=constant_stimuli,method=method, dt=dt,
                          basics_scales=basics_scales, transient=transient, seed=seed)
     
    
    pCon_reg = 0.1
    pulse0=(1100, 1200)
    pulse1=(2100, 2200)
    rate=30
    gmax_reg = 0.1
    pfail_reg=0
    PC_L23 = cortex.neuron_idcs(('PC_L23',0))
    PC_L5 = cortex.neuron_idcs(('PC_L5',0))
    
    cortex.set_regular_stimuli('poisson1', 100, ['AMPA', 'NMDA'], PC_L23, pCon=pCon_reg, rate=rate, start=pulse0[0], stop=pulse0[1], gmax=gmax_reg, pfail=pfail_reg)
    cortex.set_regular_stimuli('poisson1', 100, ['AMPA', 'NMDA'], PC_L5, pCon=pCon_reg, rate=rate, start=pulse1[0], stop=pulse1[1], gmax=gmax_reg, pfail=pfail_reg)

    cortex.run(duration)
    
    Figures.fig17(cortex, pulse0, pulse1, simulation_dir)

@time_report('Task set')
def TaskSet(simulation_dir, seed=None):
    task1(simulation_dir=simulation_dir, seed=seed)
    print(('='*50 + '\n')*2)
    task2(simulation_dir)
    print(('='*50 + '\n')*2)
    task3(simulation_dir=simulation_dir, seed=seed)
    print(('='*50 + '\n')*2)
    task4(simulation_dir=simulation_dir, seed=seed)
    print(('='*50 + '\n')*2)
    task5(simulation_dir=simulation_dir, seed=seed)
    print(('='*50 + '\n')*2)
    task6(simulation_dir=simulation_dir, seed=seed)
    print(('='*50 + '\n')*2)
    task7(simulation_dir=simulation_dir, seed=seed)
    print(('='*50 + '\n')*2)
    task8(simulation_dir=simulation_dir, seed=seed)
    print(('='*50 + '\n')*2)
    

if __name__=='__main__':
    
    simulation_dir = set_simulation_dir()
    seed = 0

    # task1(simulation_dir=simulation_dir, seed=seed)
    # task2(simulation_dir)
    # task3(simulation_dir=simulation_dir, seed=seed
    # task4(simulation_dir=simulation_dir, seed=seed)
    # task5(simulation_dir=simulation_dir, seed=seed)
    # task6(simulation_dir=simulation_dir, seed=seed)
    # task7(simulation_dir=simulation_dir, seed=seed)
    # task8(simulation_dir=simulation_dir, seed=seed)
    TaskSet(simulation_dir=simulation_dir, seed=seed)
