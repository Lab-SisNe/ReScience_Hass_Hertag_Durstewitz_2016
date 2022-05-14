import numpy as np
from scipy.signal import argrelextrema as argex
from scipy.signal import periodogram
import os
from scipy.ndimage import gaussian_filter1d as gf1d
from SimulationSetup import *

simulation_dir = set_dir()

time = np.genfromtxt('Original_t.txt', delimiter=';')
Itot = np.genfromtxt('Original_I.txt', delimiter=';')



time_step = 0.05

fs = 1000/time_step

t_arr = time[time>=1000]-1000
I_arr = Itot[time>=1000]

t_arr_bin = (t_arr//time_step).astype(int)



I_arr_final = np.zeros(np.max(t_arr_bin)+1)
t_arr_final = np.arange(0, (np.max(t_arr_bin)+1)*time_step, time_step)


for i in range(np.max(t_arr_bin)+1):
    
    I_arr_final[i] = np.mean(I_arr[t_arr_bin==i])
 
freq_arr, power_arr = periodogram(I_arr_final, fs)



log_power = np.log(power_arr)
log_freq = np.log(freq_arr)

filtered_log_power = gf1d(log_power, 11)
filtered_log_freq = log_freq



fig, ax = subplots(figsize=(12,10))
plt.xlim(2, 7)
ax.set_xlabel('log(Frequency) (log[Hz])', fontsize=26)
ax.set_ylabel('log(Power) (arbitrary unit)', fontsize=26)
ax.plot(filtered_log_freq, filtered_log_power)
ax.xaxis.set_tick_params(labelsize=26)
ax.yaxis.set_tick_params(labelsize=26)
plt.savefig('{}/Fqspectrum_filtered.png'.format(simulation_dir))
