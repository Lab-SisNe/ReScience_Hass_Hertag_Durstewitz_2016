#################################################################################################################
Simulation codes from Rempel, Kamiji and Roque (currently under review), "[RE] A Detailed Data-Driven Network Mo-
del of Prefrontal cortex Reproduces Key Features of In Vivo Activity", developed as a replication of Hass, Hertaeg
and Durstewitz (2016) [1].
#################################################################################################################

This set of files enables one to run and analyse the simulation results of the reimplemented prefrontal cortex mo-
del.

################################## FILES ##################################
To perform a simulation, one should run one of the 'Main' files. 'Main.py' is the default main script, where details
of simulation and analyses are specified, as described in the next section. Other main files are set to reproduce
specific analyses:

- 'Main_paramdistr.py': distribution of membrane parameters throughout the cell types in the column.
- 'Main_spikingneurons.py': counting of spiking neurons and comparison of connectivity and membrane and synaptic pa-
rameters between spiking and not spiking neurons.
- 'Main_stanalysis.py': spike-trains analyses (ISI mean and CV, correlations) and V analyses.
- 'Main_LFP.py': local field potential analyses.
- 'Main_stanalysisLFP.py': this script performs the same as 'Main_stanalysis.py' and 'Main_LFP.py'.
- 'Main_regularstimulation.py': regular stimulation of the network with different values for std of membrane parame-
ters.
- 'Main_poissonstimulation.py': Poisson stimulation of the network with original and 40% of original inhibitory
g_syn values.
- 'Main_LFP_rk2.py': the same as 'Main_LFP.py' using Brian 2 'gsl_rk2' algorithm.
- 'Main_stanalysis_rk2.py': the same as 'Main_stananalysis.py' using Brian 2 'gsl_rk2' algorithm.
- 'Main_stanalysisLFP_rk2.py': the same as 'Main_stanalysisLFP.py' using Brian 2 'gsl_rk2' algorithm.
 
'ParamSetup.py' sets data structures of general network features, as the relative size of neuron groups, connection
 probabilities and distribution of membrane and synaptic parameters. 'NetworkSetup.py' generates connection matrix
and parameter sets for the specific network instance to be run. 'CortexSetup.py' initializes Brian 2 objects with 
the structure previously obtained. 'SimulationSetup.py' complements the simulation setup. 'AuxiliarEquations.py' con-
tains equations that are necessary to calculate f-I functions and to draw phase portraits. 'SingleNeuron.py' allows
one to simulate single simpAdEx neurons.

################################## SIMULATION SETUP ##################################
In each run, a simulation directory is created as 'Simulation_{}', where {} is the lowest not used integer. The files
thar are automatically created during the simulation are saved to the simulation directory.
At the beginning of simulation, details of simulation setup are saved as 'SIMULATION_INFO.txt'.
Parameters and connections setup can be generated and saved as txt files for later use.  

#################### Control parameters ####################
* 'control_params' sets the main simulation control features:
- Duration: duration of simulation in ms.
- Time step: step for numerical integration in ms.
- Method: numerical integration method (as in Brian 2, i.e. 'rk4', 'gsl_rk2').
- Neurons per stripe: size of each column (to be rounded due to subgroups fracions)
- Stripes: number of columns
- Recover/save: if 'save', the network connectivity and parameters are generated and saved. If X (integer), network se-
tup is retrieved from 'Simulation_X' rather than generated from scratch. If False, network connectivity and parameters
 are generated and not saved.
- Run: if True, the simulation is performed after network setup. Otherwise simulation is not performed (this mode can
be used to generate and save network setups for later use, or else to analyse connectivity and parameter distributions).
- Seed: if Y (float/integer), a seed for random number generator only in the proper simulation with Brian 2 (not in the
 network setup) is set.

#################### Scales ####################
g_max_scale, pCon_scale and param_std_scale deal with network parameter setup and are saved with it. In order to change
these scales in a new simulation, they must not be recovered from a previous one.

* g_max_simulation: scale factor for g_max for (in order):
- excitatory to excitatory synapses
- inhibitory to excitatory synapses
- excitatory to inhibitory synapses
- inhibitory to inhibitory synapses
Default: [1, 1, 1, 1]

-----------------------------------------------

* pCon_scale: scale for pCon (in order):
- excitatory to excitatory synapses
- inhibitory to excitatory synapses
- excitatory to inhibitory synapses
- inhibitory to inhibitory synapses
Default: [1, 1, 1, 1]

-----------------------------------------------

* param_std_scale: scale for membrane parameters std in each group
10 x 14 matrix
Order in axis 0: [0] C // [1] g_L // [2] E_L // [3] delta_T // [4] V_up //
                 [5] tau_w //[6] a // [7] b // [8] V_r // [9] V_T // 

Order in axis 1: [0] PC_L23 // [1] IN_L_L23 // [2] IN_L_d_L23 // [3] IN_CL_L23 // [4] IN_CL_AC_L23 // [5] IN_CC_L23 // [6] IN_F_L23 //
                  [7] PC_L5 // [8] IN_L_L5 // [9] IN_L_d_L5 // [10] IN_CL_L5 // [11] IN_CL_AC_L5 // [12] IN_CC_L5 // [13] IN_F_L5 //

Default: 10 x 14 matrix full of ones

-----------------------------------------------

[source/target][AMPA/GABA/NMDA]_gmaxscale work inside CortexNetwork, and changes do take place if network setup is retrieved
from a previous one.
* source[AMPA/GABA/NMDA]_gmaxscale: scale synapses where the specified group is pre-synaptic.
* target[AMPA/GABA/NMDA]_gmaxscale: scale synapses where the specified group is post-synaptic.

Structure: [source/target][AMPA/GABA/NMDA]_gmaxscale = list of lists
Outer list: each scale to implement
Inner list: scale details: [group, stripe, scale value]

Ex. targetAMPA_gmaxscale = [['all', 0, 0.2]]

#################### Clustering ####################

* recur_clustering: if neighborhood rule is to be applied or not

#################### Stimuli ####################
If no specifications, default settings are applied

* constant_current: corresponds to the background current (in pA)
List of lists:
Outer list: list of configuration for each column
Inner list: configuration for each group inside the column

Default: [[250, 200, 200, 200, 200, 200, 200, 250, 200, 200, 200, 200, 200, 200]]

-----------------------------------------------

* fluctuating_current: further fluctuating current
Default: no fluctuating current

Example
fluctuating_current = [0, # Start
                        1000, # End
                        [['sin(t)' for i in range(14)] for k in range(control_param['Stripes'])] # one function to each group and stripe, written as str, as a function on t
                        ]

-----------------------------------------------

* PoissonStimuli: Poisson stimulation
Outermost structure: List of lists
Main list: all sets of stimuli
Sublist: Individual set of stimulus
In each sublist:   [0] - Number of cells in the poissonic spiking group (CV ~ 1)
                   [1] - Excitatory(1) or inhibitory (2)
                   [2] - Frequence (in Hz)
                   [3] - Synaptic strength (in nS)
                   [4] - failure probability
                   [5] - Start time (in ms)
                   [6] - End time (in ms)
                   [7] - Innermost list of lists: targets
                         main list: set of targets
                         sublist: each target ([0] - group index; [1] - stripe index; [2] - connection prob)

Example:
PoissonStimuli = [
                    [100, 1, 30, 2, 0,  1000, 1100, [[0, 0, 0.1],],],
                    [100, 1, 30, 2, 0,  2500, 2600, [[7, 0, 0.1],],],
                  ]

-----------------------------------------------

* Regular stimuli
Outermost structure: List of lists
Main list: all sets of stimuli
Sublist: Individual set of stimulus
In each sublist:   [0] - Number of cells in the regular spiking group (CV = 0)
                   [1] - Excitatory(1) or inhibitory (2)
                   [2] - Number of spikes
                   [3] - Synaptic strength (in nS)
                   [4] - failure probability
                   [5] - Start time (in ms)
                   [6] - End time (in ms)
                   [7] - Innermost list of lists: targets
                         main list: set of targets
                         sublist: each target ([0] - group index; [1] - stripe index; [2] - connection prob)
                        
Example:
RegularStimuli = [
                 [1, 1, 250, 0.1, 0, 1000, 1005, [[0, 0, 0.1],],],                
                 [1, 1, 500, 0.1, 0, 2500, 2505, [[0, 0, 0.1],],],
                ]

#################### Monitors ####################
* NeuronMonitors: monitors to record neuron variables
Outermost structure: List of lists
Main list: all sets of neuron monitors
Sublist: Individual set of neuron monitor
In each sublist: [0] - Dictionary key to the monitor
                 [1] - List of recording variables
                 [2] - Innermost list of list
                       Main list: all recording neuron groups
                       Sublist: Each recording neuron group ([0] - Group name; [1] - Stripe index)
                       (If [['all', 'all'],]: all neurons are recorder)

   
Example:
NeuronMonitor = [
                ['V', ['V'], [['all', 'all'],],], #--> to V analysis
                ['w', ['w'], [['all', 'all'],],],             
                ['I_tot', ['I_tot',], [['all', 'all'],],], #--> to LFP
                ['I_AMPA', ['I_AMPA',], [['all', 'all'],],],
                ['I_NMDA', ['I_NMDA',], [['all', 'all'],],],
                ['I_GABA', ['I_GABA',], [['all', 'all'],],],
                  ]

-----------------------------------------------

* SynapsesMonitors: monitors to record synapse variables
Outermost structure: List of lists
Main list: all sets synapse monitors
Sublist: Individual set of synapse monitor
In each sublist: [0] - Dictionary key to the monitor
                 [1] - List of recording variables
                 [2] - Innermost list of list
                       Main list: all recording synapse groups
                       Sublist: Each recording synapse group (defined by the target and the source neuron groups)
                       ([0] - Target group index; [1] - Target stripe index; [2] - Source group index; [3] - Source stripe index)

Example:
SynapsesMonitor = [
                  ['R', ['R'], [['PC_L23', 0, 'PC_L5', 0],],],        
                  ]   


################################## ANALYSES SETUP ##################################

analysis_params is a dictionary containing automatic analyses specifications in each of its items.
Outputs are pictures and reports saved in subdirectories inside 'Simulation_X' and lists with results for further manipulation if desired.

* 'ISI' - analyses of spike trains - zero-lag cross-correlations, autocorrelations (calculated according to formula extracted from [2] and ISI statistics
'Group': groups to analyse - (list of lists)
                             Main list: set of groups
                             Sublist: group ([0] - Group name; [1] - stripe index)
'Start': start time (in ms)
'Stop': stop time (in ms)
'Time bin': size of Time bin (in ms)
'Minimum spike number': minimum number of spikes in the selected neurons

Output:
spiketime_list elements: list of spike times
zerolagCC_list elements: list of the pairwise zero-lag crosscorrelation between the analysed neurons
ISImean_list elements: list of the ISI means for each analysed neuronc
ISIstd_list elements: list of the ISI std for each analysed neuron
ISICV_list elements: list of the ISI std for each analysed neuron
binned_spiketimes_list elements: 2d-array with time bins in each row (0: bin with no spike; 1: bin with spike)
return_expecvalue_list elements: list of binned spiking expected value for each neuron (# spikes / # bins)
return_var_list elements: list of binned spiking variance for each neuron
autocorrvalue_list elements: 2d array with autocorrelation for each neuron in the rows
autocorrt_list elements: list of time values used in autocorrelation

-----------------------------------------------

* 'DAcorrelation' - correlations using formula extracted from [3]

'Group': groups to analyse - (list of lists)
                             Main list: set of groups
                             Sublist: group ([0] - Group name; [1] - stripe index)
'Start': start time (in ms)
'Stop': stop time (in ms)
'Time bin': size of time bins (in ms)
'Minimum spike number': minimum number of spikes in the selected neurons

Output:
zerolagCC_list elements: list of the pairwise zero-lag crosscorrelation between the analysed neurons
autocorrvalue_list elements: 2d array with autocorrelation for each neuron in the rows
autocorrt_list elements: list of time values used in autocorrelation

-----------------------------------------------

* 'V' - analyses of V traces after extraction of spike events (individual mean, std and V_T - mean)

'Group': groups to analyse - (list of lists)
                             Main list: set of groups
                             Sublist: group ([0] - Group name; [1] - stripe index)
'Start': start time (in ms)
'Stop': stop time (in ms)
'Minimum spike number': minimum number of spikes in the selected neurons

Output:
Vmean_list elements: list of individual V mean
Vstd_list elements: list of individual V std
Vsubthres_list elements: list of individual V_T - V mean

-----------------------------------------------

* 'Vcorr' - analysis of V traces without spike event extraction and V correlations

'Group': groups to analyse - (list of lists)
                             Main list: set of groups
                             Sublist: group ([0] - Group name; [1] - stripe index)
'Start': start time (in ms)
'Stop': stop time (in ms)
'Minimum spike number': minimum number of spikes in the selected neurons

Output:
monitor_t: list of time points
Vindividualstd_list elements: lists of std of individual V traces
Vindividualmean_list elements: lists of mean of individual V traces
Vgroup_list elements: collective V traces (as mean of inidividual traces at each time point) 
VzerolagCC_list elements: lists of pairwise zero-lag cross-correlation between V traces (without normalization)
VnormalizedzerolagCC_list: lists of pairwise zero-lag cross-correlation between V traces (normalized by stds )

-----------------------------------------------

* 'Frequency' - analysis of SPD of currents. Total I_tot in the column is used as estimation of LFP.

'Group': groups to analyse - (list of lists)
                             Main list: set of groups
                             Sublist: group ([0] - Group name; [1] - stripe index)
'Source': current component to analyse (i.e. I_tot; I_GABA; I_NMDA; I_AMPA; I_EXC; I_inj)
'Start': start time (in ms)
'Stop': stop time (in ms)
'Minimum spike number': minimum number of spikes in the selected neurons
'Maximum spike number': maximum number of spikes in the selected neurons

Output
Imonitort_list elements: lists of time points 
I_list elements: lists of I values
LFPfrequency_list: lists of frequency values
LFPpower_list: lists of the corresponding power values

-----------------------------------------------

* 'Populational rate' - Rate as the number of spikes in the group by the number of neurons in the spike and by the time interval in each bin

'Group': groups to analyse - (list of lists)
                             Main list: set of groups
                             Sublist: group ([0] - Group name; [1] - stripe index)
'Source': current component to analyse (i.e. I_tot; I_GABA; I_NMDA;
'Time bins': size of time bins (in ms)

Output
Figure Pop_rate_{}.png

popratet_lists elements: lists of time points 
popratecount_lists elements: lists of the corresponding spike counts
popratefreq_lists elements: lists of the corresponding spike counts

-----------------------------------------------

* 'Rate distribution'- 	Proportion of neurons spiking in each frequency band (defined in 'Bins')

'Group': groups to analyse - (list of lists)
                             Main list: set of groups
                             Sublist: group ([0] - Group name; [1] - stripe index)
'Start': start time (in ms)
'Stop': stop time (in ms)
'Bins': list of partitioning values of the frequency bands

Output
File: Report_{}_{}.txt
ratedistribution_total_list elements: total number of cells
ratedistribution_count_list elements: number of cells in each frequency band
ratedistribution_neuron_list: neuron indices in each frequency band

-----------------------------------------------

Raster plots can be automatically built and saved with raster_plot method of CortexNetwork.

Example:
cortex.raster_plot(simulation_dir, tlims=[max(0, control_param['Duration']-4000), control_param['Duration']])
    

-----------------------------------------------
Reference:
[1] J. Hass, L. Hertäg, and D. Durstewitz. “A Detailed Data-Driven Network Model of Prefrontal Cortex Reproduces
Key Features of In Vivo Activity.” In: PLOS Computational Biology 12.5 (May 2016), pp. 1–29. DOI: 10.1371/jour-
nal.pcbi.1004930. URL: https://doi.org/10.1371/journal.pcbi.1004930.

[2] C. S. Quiroga-Lombard, J. Hass, and D. Durstewitz. “Method for stationarity-segmentation of spike train
data with application to the Pearson cross-correlation.” In: Journal of Neurophysiology 110.2 (2013). PMID:
23636729, pp. 562–572. DOI: 10.1152/jn.00186.2013. 

[3] P. Dayan and L. F. Abbott. Theoretical Neuroscience: Computational and Mathematical Modeling of Neural
Systems. Cambridge, MA: The MIT Press, 2005.

----------------------------------------------------------------------------------------------

For questions of comments, please contact:
Marcelo Rafael Silva Rempel (marcelorempel@usp.br / marcelo.rafael@hotmail.com)

Department of Physics
Faculty of Philosophy, Sciences and Letters at Ribeirão Preto
University of Sao Paulo
