import numpy as np
import os
from matplotlib import pyplot as plt
from PearsonCorrelation import *
from SimulationSetup import *

simulation_dir = set_dir()

stset = []
Stfile = 'spikes.txt'

with open('Original_spiketime.txt', 'r') as f:
    tex = f.read()

spikingneurons = tex.split(';')
spike_list = []
for i in range(len(spikingneurons)):
    _ = np.asarray(spikingneurons[i].split(',')).astype(float)
    spike_list.append(_)

for i in range(len(spike_list)):
    _ = spike_list[i]
    _= _[_>=1000]
    if len(_)>20:
        stset.append(_)
        
ISI_mean = []
ISI_std = []
ISI_CV = []                

for i in range(len(stset)):
   
    T_spike = stset[i]
    ISI = np.diff(T_spike)
    ISI_mean.append(np.mean(ISI))
    ISI_CV.append(np.std(ISI)/np.mean(ISI))
    ISI_std.append(np.std(ISI))

ISI_file = '{}/stanalysis.txt'.format(simulation_dir)

with open(ISI_file, 'a') as f:
     
    print('Spiking neurons (n spikes >= {}):'.format(21), len(spike_list), end='\n\n', file=f)
    
    print('ISI mean mean:', np.mean(ISI_mean), file=f)
    print('ISI mean std:', np.std(ISI_mean), end='\n\n', file=f)
    
    print('ISI std mean:', np.mean(ISI_std), file=f)
    print('ISI std std:', np.std(ISI_std), end='\n\n', file=f)
    
    print('ISI CV mean:', np.mean(ISI_CV), file=f)
    print('ISI CV std:', np.std(ISI_CV), end='\n\n', file=f)

fig, [ax0, ax1, ax2] = plt.subplots(1,3, figsize=(18, 10))


fig.text(0.01, 0.5, 'relative frequence', va='center', rotation='vertical', fontsize=26)
fig.text(0.315, 0.5, 'relative frequence', va='center', rotation='vertical', fontsize=26)
fig.text(0.62, 0.5, 'cross-correlation', va='center', rotation='vertical', fontsize=26)

fig.subplots_adjust(left=0.08,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.35, 
                    hspace=0.1)

lag, correlation = Pearson_correlation_group(stset, 1000, 61000, 2, -30, 30)
    

ax0.hist(ISI_mean, bins=int(np.sqrt(len(ISI_mean))), color='blue')
ax0.set_xlim(-25, 2250)
plt.sca(ax0)
plt.xticks([0, 1000, 2000], [0, 1, 2], fontsize=26)
yt1 = np.asarray([0.25, 0.50])
yt0 = yt1 * len(ISI_mean)               
plt.yticks(yt0, yt1)
yt = np.asarray([50.5, 101])
plt.yticks(yt0, yt1, fontsize=26)
ax0.set_xlabel('ISI mean (s)', fontsize=26)

ax1.hist(ISI_CV, bins=int(np.sqrt(len(ISI_CV))), color='blue')
plt.sca(ax1)
ax1.set_xlim(-0.2, 4.2)
plt.xticks([0, 1, 2, 3, 4], fontsize=26)
ax1.set_xlabel('ISI CV', fontsize=26)
yt1 = np.asarray([0.1, 0.2])
yt0 = yt1 * len(ISI_CV)               
plt.yticks(yt0, yt1, fontsize=26)

lag = lag/1000
ax2.plot(lag, correlation, color='blue')
ax2.set_xlim(-0.030, 0.030)
ax2.set_ylim(-0.0004, 0.0012)
plt.sca(ax2)
plt.xticks([-0.015, 0, 0.015], [-0.015, '0', 0.015], fontsize=26)
yt1 = np.asarray([0, 0.001])
           
plt.yticks(yt1,['0', '0.001'], fontsize=26)
ax2.set_xlabel('lag (s)', fontsize=26)

fig.savefig('{}/spanalysis.png'.format(simulation_dir))
    