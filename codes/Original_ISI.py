import numpy as np
import os
from matplotlib import pyplot as plt
from AuxiliarEquations import *
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
        
ISImean = []
ISIstd = []
ISICV = []                

for i in range(len(stset)):
   
    T_spike = stset[i]
    ISI = np.diff(T_spike)
    ISImean.append(np.mean(ISI))
    ISICV.append(np.std(ISI)/np.mean(ISI))
    ISIstd.append(np.std(ISI))

ISI_file = '{}/stanalysis.txt'.format(simulation_dir)

with open(ISI_file, 'a') as f:
     
    print('Spiking neurons (n spikes >= {}):'.format(21), len(spike_list), end='\n\n', file=f)
    
    print('ISI mean mean:', np.mean(ISImean), file=f)
    print('ISI mean std:', np.std(ISImean), end='\n\n', file=f)
    
    print('ISI std mean:', np.mean(ISIstd), file=f)
    print('ISI std std:', np.std(ISIstd), end='\n\n', file=f)
    
    print('ISI CV mean:', np.mean(ISICV), file=f)
    print('ISI CV std:', np.std(ISICV), end='\n\n', file=f)

lag, correlation = Pearson_correlation_group(stset, 1000, 61000, 2, -30, 30)


fig, [ax0, ax1, ax2] = plt.subplots(1,3, figsize=(18, 10))

fig.text(0.01, 0.5, 'relative frequency', va='center', rotation='vertical', fontsize=26)
fig.text(0.315, 0.5, 'relative frequency', va='center', rotation='vertical', fontsize=26)
fig.text(0.62, 0.5, 'cross-correlation', va='center', rotation='vertical', fontsize=26)

fig.subplots_adjust(left=0.08,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.35, 
                    hspace=0.1)


ax0.hist(ISImean, bins=int(np.sqrt(len(ISImean))), color='blue')
ax0.set_xlim(-25, 2250)
plt.sca(ax0)
plt.xticks([0, 1000, 2000], [0, 1, 2], fontsize=26)
yt1 = np.asarray([0.25, 0.50])
yt0 = yt1 * len(ISImean)               
plt.yticks(yt0, yt1)
yt = np.asarray([50.5, 101])
plt.yticks(yt0, yt1, fontsize=26)
ax0.set_xlabel('ISI mean (s)', fontsize=26)
x0, x1 = ax0.get_xlim()
y0, y1 = ax0.get_ylim()
ax0.text(0.15*x0 + 0.85*x1, 0.05*y0 + 0.95*y1, '(a)', fontsize=26)


ax1.hist(ISICV, bins=int(np.sqrt(len(ISICV))), color='blue')
plt.sca(ax1)
ax1.set_xlim(-0.2, 4.2)
plt.xticks([0, 1, 2, 3, 4], fontsize=26)
ax1.set_xlabel('ISI CV', fontsize=26)
yt1 = np.asarray([0.1, 0.2])
yt0 = yt1 * len(ISICV)               
plt.yticks(yt0, yt1, fontsize=26)
x0, x1 = ax1.get_xlim()
y0, y1 = ax1.get_ylim()
ax1.text(0.15*x0 + 0.85*x1, 0.05*y0 + 0.95*y1, '(b)', fontsize=26)

lag = lag/1000
ax2.plot(lag, correlation, color='blue')
ax2.set_xlim(-0.030, 0.030)
ax2.set_ylim(-0.0002, 0.0012)
plt.sca(ax2)
plt.xticks([-0.015, 0, 0.015], [-0.015, '0', 0.015], fontsize=26)
yt1 = np.asarray([0, 0.001])
plt.yticks(yt1,['0', '0.001'], fontsize=26)
ax2.set_xlabel('lag (s)', fontsize=26)
x0, x1 = ax2.get_xlim()
y0, y1 = ax2.get_ylim()
ax2.text(0.15*x0 + 0.85*x1, 0.05*y0 + 0.95*y1, '(c)', fontsize=26)

if not os.path.isdir('Figures'):
    os.mkdir('Figures')
fig.savefig('Figures/Fig06.png'.format(simulation_dir))
 
CC0 = Zero_lag_Pearson_correlation_group(stset, 1000, 61000, 2)
meanCC0 = np.mean(CC0)
stdCC0 = np.std(CC0)

with open('{}/zero_lag_Pearson_CC.txt'.format(simulation_dir), 'a') as f:
    print('Zero lag Pearson cross-correlation\n', file=f)
    print('Mean:', meanCC0, file=f)
    print('SD:', stdCC0, file=f)


time_bin=2
Nneurons = len(stset)
time_interval = 60000
N_bins = int(round(time_interval/time_bin, 0))
binned_spiketimes = np.zeros((Nneurons, N_bins))

auto_bins = N_bins//2

auto_cov_arr = np.zeros((Nneurons, auto_bins))
expecvalue_list = []
var_list = []
r_list = []
sp_arr = []

for i in range(len(stset)):
    T_spike = np.asarray(stset[i]) - 1000
    
    T1 = np.floor(T_spike/time_bin).astype(int)

    count = len(T1)                        
    expecvalue = count/N_bins
    expecvalue_list.append(expecvalue)
    var_list.append(expecvalue*(1-expecvalue))
    binned_spiketimes[i, T1] = 1
    sp_arr.append(T_spike)
    
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

auto_t = np.arange(auto_bins)*time_bin
meanautoC = np.mean(auto_cov_arr, axis=0)

autocorr = meanautoC[:len(auto_t)//2]
autot = auto_t[:len(auto_t)//2]

fig, [ax0, ax1] = plt.subplots(1, 2, figsize=(18,10))
fig.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.35, 
                    hspace=0.1)
ax0.plot(autot, autocorr)
ax0.set_xlim(-10, 300)
ax0.set_xlabel('lag (ms)', fontsize=26)
ax0.set_ylabel('autocorrelation', fontsize=26)
ax0.tick_params(labelsize=26)


ax1.plot(autot[1:], autocorr[1:])
ax1.set_xlim(2, 300)
ax1.set_xlabel('lag (ms)', fontsize=26)
ax1.set_ylabel('autocorrelation', fontsize=26)
plt.sca(ax1)
plt.xticks([2, 100, 200, 300], fontsize=26)
ax1.tick_params(labelsize=26)

fig.savefig('Figures/Fig07.png'.format(simulation_dir))