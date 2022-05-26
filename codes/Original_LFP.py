import numpy as np
from scipy.signal import argrelextrema as argex
from scipy.signal import periodogram
import os
from scipy.ndimage import gaussian_filter1d as gf1d
from SimulationSetup import *
from matplotlib import pyplot as plt

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


fig, ax = plt.subplots(figsize=(12,10))
ax.plot(filtered_log_freq, filtered_log_power)
ax.set_xlim(2.5, 7)
ax.set_ylim(0, 16)
ax.vlines(np.log(60), 8, 15, linestyle='--', color='black')
ax.plot([2, np.log(60)], [13 + (np.log(60)-2), 13], linestyle='--', color='blue')
ax.plot([np.log(60), 7], [12.5, 12.5 - 2*(7-np.log(60))], linestyle='--', color='blue')
ax.plot([np.log(60), 7], [10, 10 - 3*(7-np.log(60))], linestyle='--', color='blue')
ax.set_xlabel('log(Frequency) (log[Hz])', fontsize=26)
ax.set_ylabel('log(Power) (arbitrary unit)', fontsize=26)
plt.gca()
plt.yticks([0, 4, 8, 12, 16], fontsize=26)
plt.xticks([3, 4, 5, 6, 7], fontsize=26)
ax.text( 3, 7, '60 Hz', fontsize=26)
ax.arrow(3.6, 7.5, 0.4, 1, head_width=0.1)
ax.text(3.25, 14.5, '1/f', fontsize=26)
ax.text(5.5, 11, '$1/f^2$', fontsize=26)
ax.text(5.5, 3, '$1/f^3$', fontsize=26)


plt.savefig('{}/Fqspectrum_filtered.png'.format(simulation_dir))
