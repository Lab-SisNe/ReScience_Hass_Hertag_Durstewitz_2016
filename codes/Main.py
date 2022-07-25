from CortexNetwork import *
from SimulationSetup import *
from NetworkSetup import NetworkSetup
from time import time   
from warnings import filterwarnings


##### Setting simulation directory (as 'simulation_{}', with {} being the lowest integer corresponding to non-existing directory)

filterwarnings("ignore", category=RuntimeWarning) 
BrianLogger.suppress_name('resolution_conflict', 'device')

simulation_dir = set_dir()

##############################################################################
##########################                       #############################
########################## Simulation parameters #############################
##########################                       #############################
##############################################################################

################==============================###################
################===== General parameters =====###################
################==============================###################

##### *** These parameters can be changed but can't be commented out ***

#####----------||| Simulation control ||||----------#####

control_param = {'Duration': 7000, # in ms
                 'Time step': 0.05, # in ms
                 'Noise': False, # D if noise present (int or float); anything else otherwise
                 'Method': 'rk4', # brian2 integration methods
                 'Neurons per stripe': 1000,
                 'Stripes': 1,
                 'Recover/Save': False, ## For recovering: insert the directory number; for saving: insert 'Save'; else: insert False
                 'run': True, ## Insert False to avoid running; otherwise, insert True
                 'seed': None, ## Insert seed number; otherwise, insert None
                 }

###############----------||| Scales |||----------###############
##### Obs: changing in scales1 does not take place if membrane and synaptic parameters are being recovered from previous simulation
##### Order of scale factors: [0] E -> E // [1] I -> E // [2] E -> I // [3] I -> I 
gmax_scale = [1, 1, 1, 1]
pCon_scale = [1, 1, 1, 1]

##### Order in axis 0: [0] C // [1] g_L // [2] E_L // [3] delta_T // [4] V_up //
#####                  [5] tau_w //[6] a // [7] b // [8] V_r // [9] V_T // 

##### Order in axis 1: [0] PC_L23 // [1] IN_L_L23 // [2] IN_L_d_L23 // [3] IN_CL_L23 // [4] IN_CL_AC_L23 // [5] IN_CC_L23 // [6] IN_F_L23 //
#####                  [7] PC_L5 // [8] IN_L_L5 // [9] IN_L_d_L5 // [10] IN_CL_L5 // [11] IN_CL_AC_L5 // [12] IN_CC_L5 // [13] IN_F_L5 //

param_std_scale = np.ones((10, 14))# Scale of std of membrane parameters for each group

scales1 = gmax_scale, pCon_scale, param_std_scale

###### Obs: changing in scales2 can be performed in simulations with recovering of previous setup


sourceAMPA_gmaxscale = [
                        # ['all', 0, 0.3]
                        ]
sourceGABA_gmaxscale = [
                        # ['all', 0, 0.3]
                        ]
sourceNMDA_gmaxscale = [
                        # ['all', 0, 0.3]
                        ]


targetAMPA_gmaxscale = [
                        # ['IN_L', 0, 0.3],
                        
                        ]
targetGABA_gmaxscale = []
targetNMDA_gmaxscale = [
                        
                        ]

##### List of lists
##### Each scale to implement <--- scale details: [0] - target/source group; [1] - stripe; [2] - scale value
##### source AMPA/GABA/NMDA: scaling synapses from the specified source group
##### target AMPA/GABA/NMDA: scaling synapses to the specified target group

scales2 = sourceAMPA_gmaxscale, sourceGABA_gmaxscale, sourceNMDA_gmaxscale, targetAMPA_gmaxscale, targetGABA_gmaxscale, targetNMDA_gmaxscale


noise = control_param['Noise']
###############----------||| Clustering |||----------###############
##### Obs: clustering change does not take place if membrane and synaptic parameters are being recovered from previous simulation
recur_clustering = True
clustering = recur_clustering

###########=====================================================##############
###########===== Stimuli, monitors and analysis parameters =====##############
###########=====================================================##############

##### *** These parameters can be #commented out to maintain default settings ***
##### Brian units must be omitted in the parameters (i.e. constant_current[0][0] = 250, not 250*pA)
##### Current unit: pA
##### Time unit: ms


###############---------||| Stimuli |||----------###############
##### Default: constant current for one stripe (250pA in PC; 200pA in IN); no other stimuli

############### Background constant current ###############
##### List of lists (main list #Stripe  <-- sublist #Group)

Iexc = 250 
Iinh =  200 
constant_current = [[Iexc, Iinh, Iinh, Iinh, Iinh, Iinh, Iinh, Iexc, Iinh, Iinh, Iinh, Iinh, Iinh, Iinh] for i in range(control_param['Stripes'])
                    ]


############### Fluctuating current ###############

##### fluctuating_current[0]: start time (in ms)
##### fluctuating_current[1]: end time (in ms)
##### fluctuating_current[2]: list of lists (main list: #Stripe  <-- sublist: #Group)
#####                         current (in pA) as a function of t (in ms), written as string (i.e. fluctuating_current[2][0][0] = '10* sin(t*2*pi/100)')

# fluctuating_current = [0,
#                        control_param['Duration'],
#                        [['' for i in range(14)] for k in range(control_param['Stripes'])]
#                        ]

# fluctuating_current[2][0][0] = '10* sin(2*pi*t/100)'

############### Poissonic stimuli ###############

##### Outermost structure: List of lists
##### Main list: all sets of stimuli
##### Sublist: Individual set of stimulus
##### In each sublist: [0] - Number of cells in the poissonic spiking group (CV ~ 1)
####                   [1] - Excitatory(1) or inhibitory (2)
####                   [2] - Frequency (in Hz)
####                   [3] - Synaptic strength (in nS)
####                   [4] - failure probability
####                   [5] - Start time (in ms)
####                   [6] - End time (in ms)
####                   [7] - Innermost list of lists: targets
####                         main list: set of targets
####                         sublist: each target ([0] - group index; [1] - stripe index; [2] - connection prob)

PoissonStimuli = [
                    # [100, 1, 30, 2, 0,  1000, 1100, [[0, 0, 0.1],],],
                    # [100, 1, 30, 2, 0,  2500, 2600, [[7, 0, 0.1],],],
                  ]

############### Regular stimuli ###############

##### Outermost structure: List of lists
##### Main list: all sets of stimuli
##### Sublist: Individual set of stimulus
##### In each sublist: [0] - Number of cells in the regular spiking group (CV = 0)
####                   [1] - Excitatory(1) or inhibitory (2)
####                   [2] - Number of spikes
####                   [3] - Synaptic strength (in nS)
####                   [4] - failure probability
####                   [5] - Start time (in ms)
####                   [6] - End time (in ms)
####                   [7] - Innermost list of lists: targets
####                         main list: set of targets
####                         sublist: each target ([0] - group index; [1] - stripe index; [2] - connection prob)
                        

RegularStimuli = [
                # [1, 1, 250, 0.1, 0, 1000, 1005, [[0, 0, 0.1],],],                
                # [1, 1, 500, 0.1, 0, 2500, 2505, [[0, 0, 0.1],],],
                ]      


###############---------||| Monitors |||----------###############
##### Monitors others than the spiking one (which is always included) 
##### Default: none

############### Neuron Monitors ###############

##### These monitors are set in the dictionary self.neuronmonitors

##### Outermost structure: List of lists
##### Main list: all sets of neuron monitors
##### Sublist: Individual set of neuron monitor
##### In each sublist: [0] - Dictionary key to the monitor
####                   [1] - List of recording variables
####                   [2] - Innermost list of list
####                         Main list: all recording neuron groups
####                         Sublist: Each recording neuron group ([0] - Group name; [1] - Stripe index)
####                         (If [['all', 'all'],]: all neurons are recorder)


NeuronMonitor = [
                # ['V', ['V'], [['all', 'all'],],], #--> to V analysis
                # ['w', ['w'], [['all', 'all'],],],             
                # ['I_tot', ['I_tot',], [['all', 'all'],],], #--> to LFP
                # ['I_AMPA', ['I_AMPA',], [['all', 'all'],],],
                # ['I_NMDA', ['I_NMDA',], [['all', 'all'],],],
                # ['I_GABA', ['I_GABA',], [['all', 'all'],],],
                # ['g_AMPA', ['g_AMPA',], [['all', 'all'],],],
                # ['g_NMDA', ['g_NMDA',], [['all', 'all'],],],
                # ['g_GABA', ['g_GABA',], [['all', 'all'],],],
                  ]

############### Synapse Monitors ###############

##### These monitors are set in the dictionary self.synapsesmonitors

##### Outermost structure: List of lists
##### Main list: all sets synapse monitors
##### Sublist: Individual set of synapse monitor
##### In each sublist: [0] - Dictionary key to the monitor
####                   [1] - List of recording variables
####                   [2] - Innermost list of list
####                         Main list: all recording synapse groups
####                         Sublist: Each recording synapse group (defined by the target and the source neuron groups)
####                         ([0] - Target group index; [1] - Target stripe index; [2] - Source group index; [3] - Source stripe index)
####   

# SynapsesMonitor = [
#                 ['R', ['R'], [['PC_L23', 0, 'PC_L5', 0],],],        
#                  ]


###############---------||| Analysis |||----------###############

##### Automatic analysis
##### Default: none

##### Defined as a dictionary (analysis_param) with an item for each analysis type
##### Most elements of the output consists of lists;
##### each element of the list corresponds to one of the performed analysis

start = 1000
stop = control_param['Duration']


############### Interspike interval analysis ###############

##### Dicts of dicts
##### Main dict: set of ISI analysis (one for each item)
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                          Main list: set of groups
#####                                                          Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Start': start time (in ms)
#####                             'Stop': stop time (in ms)
#####                             'Time bin': size of Time bin (in ms)
#####                             'Minimum spike number': minimum number of spikes in the selected neurons

##### Output:
##### ISImean_list elements: list of the ISI means for each analysed neuronc
##### ISIstd_list elements: list of the ISI std for each analysed neuron
##### ISICV_list elements: list of the ISI std for each analysed neuron

analysis_params['ISI analysis'] = {
                            # 'Analysis 0': {'Group': [['all',0],],
                            #               'Start': start,
                            #               'Stop': stop,
                            #               'Minimum spike number': 3,
                            #               'Graphs': True,
                                          # },
                            }


############### ISI autocorrelation1 ###############

##### Dicts of dicts
##### Main dict: set of ISI autocorrelation analysis (one for each item)
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                          Main list: set of groups
#####                                                          Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Start': start time (in ms)
#####                             'Stop': stop time (in ms)
#####                             'Time bin': size of Time bin (in ms)
#####                             'Minimum spike number': minimum number of spikes in the selected neurons
#####                             'Graphs': True if graphs are to be generated; False otherwise

##### Output:
##### autocorrvalue_list elements: 2d array with autocorrelation for each neuron in the rows
##### autocorrt_list elements: list of time values used in autocorrelation

analysis_params['ISI autocorrelation1'] = {
                            # 'Analysis 0': {'Group': [['all',0],],
                            #               'Start': start,
                            #               'Stop': stop,
                            #               'Time bin': 2,
                            #               'Minimum spike number': 11,
                            #               'Graphs': False,
                            #               },
                            }

############### ISI crosscorrelation1 ###############

##### Dicts of dicts
##### Main dict: set of ISI crosscorrelation analysis (one for each item)
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                          Main list: set of groups
#####                                                          Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Start': start time (in ms)
#####                             'Stop': stop time (in ms)
#####                             'Time bin': size of Time bin (in ms)
#####                             'Minimum spike number': minimum number of spikes in the selected neurons
#####                             'Graphs': True if graphs are to be generated; False otherwise

##### Output:
##### zerolagCC_list elements: list of the pairwise zero-lag crosscorrelation between the analysed neurons

analysis_params['ISI crosscorrelation1'] = {
                            # 'Analysis 0': {'Group': [['all',0],],
                            #               'Start': start,
                            #               'Stop': stop,
                            #               'Time bin': 2,
                            #               'Minimum spike number': 11,
                            #               'Graphs': False,
                            #               },
                            }


############### ISI correlation2 #############################

##### Dicts of dicts
##### Main dict: set of ISI correlation analysis (one for each item) - with a second formula
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                          Main list: set of groups
#####                                                          Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Start': start time (in ms)
#####                             'Stop': stop time (in ms)
#####                             'Time bin': size of time bins (in ms)
#####                             'Minimum spike number': minimum number of spikes in the selected neurons
#####                             'Graphs': True if graphs are to be generated; False otherwise

##### Output:
##### DAzerolagCC_list elements: list of the pairwise zero-lag crosscorrelation between the analysed neurons
##### DAautocorrvalue_list elements: 2d array with autocorrelation for each neuron in the rows
##### DAautocorrt_list elements: list of time values used in autocorrelation

analysis_params['ISI correlation2'] = {
                                # 'Analysis 0': {'Group': [['all',0],],
                                #               'Start': start,
                                #               'Stop': stop,
                                #               'Time bin': 2,
                                #               'Minimum spike number': 11,
                                #               'Graphs': False,
                                #               },
                                    }

############### ISI Pearson correlation ###############
  
##### Dicts of dicts
##### Main dict: set of ISI Pearson cross correlation analysis (one for each item) - with a second formula
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                          Main list: set of groups
#####                                                          Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Start': start time (in ms)
#####                             'Stop': stop time (in ms)
#####                             'Time bin': size of time bins (in ms)
#####                             'Minimum spike number': minimum number of spikes in the selected neurons
#####                             'Minimum lag': minimum lag for cross-correlation analysis
#####                             'Maximum lag': maximum lag for cross-correlation analysis
#####                             'Graphs': True if graphs are to be generated; False otherwise
 
##### Output:
##### lag_array_list elements: list of correlation lags
##### correlation_mean_list elements: list of corresponding correlation values

         
analysis_params['ISI Pearson'] = {
                            # 'Analysis 0': {'Group': [['all',0],],
                            #               'Start': start,
                            #               'Stop': stop,
                            #               'Time bin': 2,
                            #               'Minimum spike number': 11,
                            #               'Minimum lag': -30,
                            #               'Maximum lag': 30,
                            #               'Graphs': True,
                                          # },
                            }



############### V analysis ###############

##### Dicts of dicts
##### Main dict: set of V analysis  (one for each item)
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                          Main list: set of groups
#####                                                          Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Start': start time (in ms)
#####                             'Stop': stop time (in ms)
#####                             'Minimum spike number': minimum number of spikes in the selected neurons
#####                             'Graphs': True if graphs are to be generated; False otherwise

##### Output:
##### Vmean_list elements: list of individual V mean
##### Vstd_list elements: list of individual V std
##### Vsubthres_list elements: list of individual V_T - V mean

analysis_params['V analysis'] = {
                        # 'Analysis 1': {'Group': [['all', 0,],],
                        #                   'Start': start,
                        #                   'Stop': stop,
                        #                   'Minimum spike number': 3,
                        #                   'Graphs': True,
                                        # },

                            }

    

############### V correlation analysis ###############

##### Dicts of dicts
##### Main dict: set of V correlation analysis (one for each item)
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                          Main list: set of groups
#####                                                          Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Start': start time (in ms)
#####                             'Stop': stop time (in ms)
#####                             'Minimum spike number': minimum number of spikes in the selected neurons
#####                             'Graphs': True if graphs are to be generated; False otherwise

##### Output:
##### monitor_t: list of time points
##### Vindividualstd_list elements: lists of std of individual V traces
##### Vindividualmean_list elements: lists of mean of individual V traces
##### Vgroup_list elements: collective V traces (as mean of inidividual traces at each time point) 
##### VzerolagCC_list elements: lists of pairwise zero-lag cross-correlation between V traces (without normalization)
##### VnormalizedzerolagCC_list: lists of pairwise zero-lag cross-correlation between V traces (normalized by stds )


analysis_params['V correlation'] = {
                        # 'Analysis 1': {'Group': [['all', 0,],],
                        #                   'Start': start,
                        #                   'Stop': stop,
                        #                   'Minimum spike number': 1,
                        #                   'Graphs': True,
                        #                 },

                            }


############### Frequency analysis ###############

##### This analysis estimates the frequency spectrum of LFP through the total synaptic current in the population
##### It's necessary to set the respective neuron monitor to the source current (i.e. neuronmonitors['I_tot']
##### recording self.group.I_tot to analyse the LFP component due to the total current)

##### Dicts of dicts
##### Main dict: set of Frequency analysis (one for each item)
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                          Main list: set of groups
#####                                                          Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Source': current component to analyse (i.e. I_tot; I_GABA; I_NMDA; I_AMPA; I_EXC; I_inj)
#####                             'Start': start time (in ms)
#####                             'Stop': stop time (in ms)
#####                             'Minimum spike number': minimum number of spikes in the selected neurons
#####                             'Maximum spike number': maximum number of spikes in the selected neurons
#####                             'Filtering': type of filtering ('Gaussian', 'Moving average'); False if none
#####                             'Filtering parameter': averaging window for moving average; sigma for gaussian filter
#####                             'Graphs': True if graphs are to be generated; False otherwise

    
##### Output
##### Imonitort_list elements: lists of time points 
##### I_list elements: lists of I values
##### LFPfrequency_list: lists of frequency values
##### LFPpower_list: lists of the corresponding power values

analysis_params['Frequency'] = {
                                # 'Analysis 1':{'Group': [['all', 0,],], #--> to LFP
                                #                         'Source': 'I_tot',
                                #                         'Start': 1000,
                                #                         'Stop': stop,
                                #                         'Minimum spike number': 0,
                                #                         'Maximum spike number': 100000,
                                #                         'Filtering': 'Gaussian',
                                #                         'Filtering parameter': 11,
                                #                         'Graphs': True,
                                #                         },                             
                                  }

                  
############### Populational rate ###############

##### Rate as the number of spikes in the group by the number of neurons in the spike and by the time interval in each bin

##### Dicts of dicts
##### Main dict: set of populational rate analysis (one for each item)
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                          Main list: set of groups
#####                                                          Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Source': current component to analyse (i.e. I_tot; I_GABA; I_NMDA; I_AMPA; I_EXC; I_inj)
#####                             'Time bins': size of time bins (in ms)



                                
analysis_params['Populational rate'] = {
                                            # 'Analysis 0': {'Group': [['PC_L23', 0],],
                                            #               'Start': start,
                                            #               'Stop': stop,
                                            #               'Time bin': 1,
                                            #               'Moving average': 31,
                                            #               'Graphs': False,
                                            #               },
                                            # 'Analysis 1': {'Group': [['PC', 0],],
                                            #               'Start': start,
                                            #               'Stop': stop,
                                            #               'Time bin': 1,
                                            #               'Moving average': 31,
                                            #               'Graphs': False,
                                            #               },

                                          }

##### Output
##### Figure Pop_rate_{}.png

##### popratet_lists elements: lists of time points 
##### popratecount_lists elements: lists of the corresponding spike counts
##### popratefreq_lists elements: lists of the corresponding spike counts

############### Rate stratification ###############

##### Proportion of neurons spiking in each frequency band (defined in 'Bins')

##### Dicts of dicts
##### Main dict: set of Frequency analysis (one for each item)
##### Subdict: Analysis setting - 'Group': groups to analyse - (list of lists)
#####                                                           Main list: set of g#####roups
#####                                                           Sublist: group ([0] - Group name; [1] - stripe index)
#####                             'Start': start time (in ms)
#####                             'Stop': stop time (in ms)
#####                             'Bins': list of partitioning values of the frequency bands

#### Output
#### File: Report_{}_{}.txt
#### ratestratification_total_list elements: total number of cells
#### ratestratification_count_list elements: number of cells in each frequency band
#### ratestratification_neuron_list: neuron indices in each frequency band
 

analysis_params['Rate stratification'] = {
                                                'Analysis 0': {'Group': [['all', 0],],
                                                                'Start': start,
                                                                'Stop':stop,
                                                                'Bins': [1/3,],
                                                                },  
                                                'Analysis 1': {'Group': [['PC', 0],],
                                                                'Start': start,
                                                                'Stop':stop,
                                                                'Bins': [1/3,],
                                                                },  
                                          
                                          }
              
################==============================###################
################===== Information saving =====###################
################==============================###################

save_info(control_param, clustering, PoissonStimuli, RegularStimuli, NeuronMonitor, 
          SynapsesMonitor, analysis_params, scales1, scales2, constant_current, 
          fluctuating_current, simulation_dir)


###################################################################
##########################            #############################
########################## Simulation #############################
##########################            #############################
###################################################################

start_time = time()

print("REPORT: Simulation directory: '{}'\n".format(simulation_dir))
print('REPORT: Starting network\n')

###############----------||| Parameters and network setup |||----------###############

if type(control_param['Recover/Save']) is int:      ##### Recover from previous simulation
    print("REPORT: Recovering membrane and synaptic parameters from 'simulation_{}'".format(control_param['Recover/Save'])) 
    NeuPar, V0, STypPar, SynPar,SPMtx, group_distr = set_setup(control_param['Recover/Save'])
    print('REPORT: Parameters recovered\n')
    
else:       ##### Generate new parameters
    print('REPORT: Generating membrane and synaptic parameters')
    NeuPar, V0, STypPar, SynPar,SPMtx, group_distr = NetworkSetup(control_param['Neurons per stripe'], control_param['Stripes'], scales1, clustering)
    print('REPORT: Parameters generated')
    if control_param['Recover/Save'] in ['Save', 'save', 'SAVE']:       #### Save parameters
        save_setup(NeuPar, V0, STypPar, SynPar,SPMtx, group_distr, simulation_dir)
        print("REPORT: Parameters saved in simulation directory ('{}')\n".format(simulation_dir))
    else:
        print('\n')

cortex = CortexNetwork(NeuPar, V0, STypPar,SynPar,SPMtx, group_distr, constant_current, fluctuating_current, scales2, noise,
                        control_param['Method'], control_param['Time step'], control_param['seed'], simulation_dir)

if type(control_param['seed']) is not bool:
    print('REPORT: Seed set to {}\n'.format(control_param['seed']))
else:   
    print('REPORT: No seed\n')

                       
minutes, seconds = divmod(time() - start_time, 60)
print("REPORT: Elapsed time for network setup: {:0>2}m {:0>2}s\n".format(int(minutes), int(round(seconds, 0))))

#####----------||| Monitors setup |||----------#####

if control_param['run']:
    if len(NeuronMonitor) + len(SynapsesMonitor): ##### Set monitors
        print('REPORT: Setting monitors')
        if len(NeuronMonitor):
            cortex.neuron_monitors(NeuronMonitor)
        if len(SynapsesMonitor):
            cortex.synapses_monitors(SynapsesMonitor) 
        print('REPORT: Monitors set\n')

###############----------||| Stimuli setup |||----------###############

if control_param['run']:
    if len(PoissonStimuli) + len(RegularStimuli)>0:
        print('REPORT: Setting stimuli protocols')
        if len(PoissonStimuli):
            cortex.poisson_input(PoissonStimuli)
        if len(RegularStimuli):
            cortex.regular_input(RegularStimuli)
        print('REPORT: Stimuli protocols set\n')

###############----------||| Actual simulation |||----------###############

if control_param['run']:
    print("REPORT: Preparing simulation")
    print('Time step: {} ms'.format(control_param['Time step']))
    print('Duration: {} ms'.format(control_param['Duration']))
    if type(control_param['Noise']) is int or type(control_param['Noise']) is float:
        print('Noise (D): {} '.format(control_param['Noise']))
    else:
        print('Noise (D): no noise')
    print('Integration method:', control_param['Method'], end='\n\n')
    
    run_time = time()
    defaultclock.dt = control_param['Duration']*ms
    seed(control_param['seed'])
    cortex.run(control_param['Duration']*ms)
     
    print()
    minutes, seconds = divmod(time() - run_time, 60)
    print("REPORT: Elapsed time for network run: {:0>2}m {:0>2}s\n".format(int(minutes), int(round(seconds, 0))))
    
    minutes, seconds = divmod(time() - start_time, 60)
    print("REPORT: Total elapsed time: {:0>2}m {:0>2}s\n".format(int(minutes), int(round(seconds, 0))))

    ###############----------||| Analysis |||----------###############
    

    if len(analysis_params['ISI analysis'].values()) + len(analysis_params['ISI autocorrelation1'].values()) + len(analysis_params['ISI crosscorrelation1'].values()) + len(analysis_params['ISI correlation2'].values()) + len(analysis_params['ISI Pearson'].values()) + len(analysis_params['V analysis'].values()) + len(analysis_params['V correlation'].values()) + len(analysis_params['Frequency'].values()) + len(analysis_params['Populational rate'].values()) + len(analysis_params['Rate stratification'].values()):
        print('REPORT: Performing analysis')
        
        if len(analysis_params['ISI analysis'].values()):
            ISImean_list, ISIstd_list, ISICV_list = cortex.ISI_analysis(analysis_params['ISI analysis'].values())           
        
        if len(analysis_params['ISI autocorrelation1'].values()):
            autocorrvalue_list, autocorrt_list = cortex.ISI_autocorrelation1(analysis_params['ISI autocorrelation1'].values())           
        
        if len(analysis_params['ISI crosscorrelation1'].values()):
            zerolagCC_list  = cortex.ISI_crosscorrelation1(analysis_params['ISI crosscorrelation1'].values())           
        
        if len(analysis_params['ISI correlation2'].values()):
            DAzerolagCC_list, DAautocorrvalue_list, DAautocorrt_list = cortex.DAcorrelation_analysis(analysis_params['ISI correlation2'].values())
        
        if len(analysis_params['ISI Pearson'].values()):
            lag_array_list, correlation_mean_list = cortex.Pearson_correlation(analysis_params['ISI Pearson'].values())
        
        if len(analysis_params['V analysis'].values()):
            Vmean_list, Vstd_list, Vsubthres_list = cortex.V_analysis(analysis_params['V analysis'].values())
        
        if len(analysis_params['V correlation'].values()):
            monitor_t,Vindividualstd_list, Vindividualmean_list, Vgroup_list, VzerolagCC_list, VnormalizedzerolagCC_list = cortex.Vcorr_analysis(analysis_params['V correlation'].values())
        
        if len(analysis_params['Frequency'].values()):
            Imonitort_list, I_list, LFPfrequency_list, LFPpower_list, filtered_LFPfrequency_list, filtered_LFPpower_list = cortex.frequency_analysis(analysis_params['Frequency'].values())
            
        if len(analysis_params['Populational rate'].values()):
            popratet_lists, popratecount_lists, popratefreq_lists, popspikescount_list = cortex.population_rate(analysis_params['Populational rate'].values())
            
        if len(analysis_params['Rate stratification'].values()):
            ratestratification_total_list, ratestratification_count_list, ratestratification_neuron_list = cortex.rate_stratification(analysis_params['Rate stratification'].values())
        
        print('REPORT: Analysis concluded\n')
    
    
    ###############----------||| Raster plot |||----------###############
    
    if not os.path.isdir('Figures'):
        os.mkdir('Figures')
    cortex.raster_plot(tlims=[max(0, control_param['Duration']-6000), control_param['Duration']], in_s=True, plotsize=5, savefig='Figures/Fig01.png')
    # cortex.raster_plot(tlims=[max(0, control_param['Duration']-1000), control_param['Duration']])
    # cortex.raster_plot(tlims=[max(0, control_param['Duration']-500), control_param['Duration']])
    
###############----------||| Parameters summary |||----------###############

GroupSetup = [[['PC_L23',0]],[['PC_L5',0]],
              [[1,0], [2,0], [8, 0], [9,0]], # IN-L and IN-L groups --> fast-spiking cells
              [[3,0], [4,0], [10,0], [11,0]], # IN-CL and IN-CL-AC groups --> bitufted cells
              [[5,0], [12,0]], # IN-CC groups --> basket cells
              [[6,0], [13,0]], # IN-F groups --> Martinotti cells
              ]
cortex.groups_params(GroupSetup)





