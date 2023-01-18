import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter1d as gf1d
import AnalysesTools as Tools
import os


def ISIfigures(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, file):
    ISImean = np.asarray(ISImean)
    ISICV = np.asarray(ISICV)
    C_lags = np.asarray(C_lags)
    crossC_mean = np.asarray(crossC_mean)
    
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
    
    plt.sca(ax0)
    plt.hist(ISImean, bins=int(np.sqrt(len(ISImean))), color='blue')
    plt.xlim(-25, 2250)
    plt.xticks([0, 1000, 2000], [0, 1, 2], fontsize=26)
    yt1 = np.asarray([0.25, 0.50])
    yt0 = yt1 * len(ISImean)               
    plt.yticks(yt0, yt1)
    yt = np.asarray([50.5, 101])
    plt.yticks(yt0, yt1, fontsize=26)
    plt.xlabel('ISI mean (s)', fontsize=26)
    x0, x1 = ax0.get_xlim()
    y0, y1 = ax0.get_ylim()
    plt.text(0.15*x0 + 0.85*x1, 0.05*y0 + 0.95*y1, '(a)', fontsize=26)
    
    plt.sca(ax1)
    plt.hist(ISICV, bins=int(np.sqrt(len(ISICV))), color='blue')
    plt.xlim(-0.2, 4.2)
    plt.xticks([0, 1, 2, 3, 4], fontsize=26)
    plt.xlabel('ISI CV', fontsize=26)
    yt1 = np.asarray([0.1, 0.2])
    yt0 = yt1 * len(ISICV)               
    plt.yticks(yt0, yt1, fontsize=26)
    x0, x1 = ax1.get_xlim()
    y0, y1 = ax1.get_ylim()
    plt.text(0.15*x0 + 0.85*x1, 0.05*y0 + 0.95*y1, '(b)', fontsize=26)
    
    plt.sca(ax2)
    C_lags0 = C_lags[C_lags>=0]
    crossC_mean0 = crossC_mean[C_lags>=0]
    
    C_lags1=C_lags0[-1:0:-1]
    C_lags2 = np.concatenate(( -C_lags1, C_lags0))
    
    crossC_mean1 = crossC_mean0[-1:0:-1]
    crossC_mean2 = np.concatenate((crossC_mean1, crossC_mean0))
    
    C_lags2 = C_lags2/1000
    
    if Correlation_sigma>0:
        ax2.plot(C_lags2, gf1d(crossC_mean2,Correlation_sigma), color='blue')
    else:
        ax2.plot(C_lags2, crossC_mean2, color='blue')
        
    ax2.set_xlim(-0.030, 0.030)
    # ax2.set_ylim(-0.0002, 0.0012)
    plt.xticks([-0.015, 0, 0.015], [-0.015, '0', 0.015], fontsize=26)
    yt1 = np.asarray([0, 0.001])
    plt.yticks(yt1,['0', '0.001'], fontsize=26)
    plt.xlabel('C_lags (s)', fontsize=26)
    x0, x1 = ax2.get_xlim()
    y0, y1 = ax2.get_ylim()
    plt.text(0.15*x0 + 0.85*x1, 0.05*y0 + 0.95*y1, '(c)', fontsize=26)
      
    
    plt.savefig(file)
    
    plt.show()
    plt.close()

def autoCfigures(C_lags, autoC_mean, file):
    fig, [ax0, ax1] = plt.subplots(1, 2, figsize=(18,10))
    fig.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.35, 
                        hspace=0.1)
    
    ax0.plot(C_lags, autoC_mean)
    ax0.set_xlim(-10, 300)
    ax0.set_xlabel('lag (ms)', fontsize=26)
    ax0.set_ylabel('autocorrelation', fontsize=26)
    ax0.tick_params(labelsize=26)
    
    
    ax1.plot(C_lags[1:], autoC_mean[1:])
    ax1.set_xlim(2, 300)
    ax1.set_xlabel('lag (ms)', fontsize=26)
    ax1.set_ylabel('autocorrelation', fontsize=26)
    plt.sca(ax1)
    plt.xticks([2, 100, 200, 300], fontsize=26)
    ax1.tick_params(labelsize=26)
    fig.savefig(file)
    
def SPDfigures(frequency, power, file):
    fig, ax = plt.subplots(figsize=(12,10))
    power = power
    ax.plot(frequency, power)
    ax.set_xlim(2.5, 7)
    ax.set_ylim(-2, 18)
    ax.vlines(np.log(60), 6, 15, linestyle='--', color='black')
    ax.plot([2, np.log(60)], [13.5 + (np.log(60)-2), 13.5], linestyle='--', color='blue')
    ax.plot([np.log(60), 7], [12.5, 12.5 - 2*(7-np.log(60))], linestyle='--', color='blue')
    ax.plot([np.log(60), 7], [8, 8 - 3*(7-np.log(60))], linestyle='--', color='blue')
    ax.set_xlabel('log(Frequency) (log[Hz])', fontsize=26)
    ax.set_ylabel('log(Power) (arbitrary unit)', fontsize=26)
    plt.gca()
    plt.yticks([0, 4, 8, 12, 16], fontsize=26)
    plt.xticks([3, 4, 5, 6, 7], fontsize=26)
    ax.text( 3, 7, '60 Hz', fontsize=26)
    ax.arrow(3.6, 7.5, 0.4, 1, head_width=0.1)
    ax.text(3.25, 15.5, '1/f', fontsize=26)
    ax.text(5.5, 11, '$1/f^2$', fontsize=26)
    ax.text(5.5, 1, '$1/f^3$', fontsize=26)

    fig.savefig(file)
    
def fig01(cortex, xlim, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    Tools.raster_plot(cortex, xlim=xlim, savefig=path+'\\Figures\\Fig01.png')

def fig02(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    ISIfigures(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, path+'\\Figures\\Fig02.png')
    
def fig03(C_lags, autoC_mean, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    autoCfigures(C_lags, autoC_mean, path+'\\Figures\\Fig03.png')
    
def fig04(frequency, power, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    SPDfigures(frequency, power, path+'\\Figures\\Fig04.png')
   
def fig05(VstdALL, VstdPC, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
     
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.set_xlabel('V standard deviation (mV)', fontsize=26)
    ax.set_ylabel('relative frequency', fontsize=26)
    h = ax.hist(VstdALL, bins=int(round(np.sqrt(len(VstdALL)), 0)))
    ax.xaxis.set_tick_params(labelsize=26)
    y1 = np.max(h[0])/len(VstdALL)
    if y1 % 0.1 >= 0.05:
        yf = 0.1*(int(y1*10) + 1)
    else:
        yf = 0.1*int(y1*10)

    yarr1 = np.arange(0, yf+0.01, 0.1)
    yarr0 = yarr1*len(VstdALL)

    plt.gca()
    plt.yticks(yarr0, np.round(yarr1, 1), fontsize=26)

    ax.yaxis.set_tick_params(labelsize=26)
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.text(0.15*x0 + 0.85*x1, 0.05*y0 + 0.95*y1, '(a)', fontsize=26)

    
    fig.savefig(path+'\\Figures\\Fig05a.png')

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.set_xlabel('V standard deviation (mV)', fontsize=26)
    ax.set_ylabel('relative frequency', fontsize=26)
    h = ax.hist(VstdPC, bins=int(round(np.sqrt(len(VstdPC)), 0)))
    ax.xaxis.set_tick_params(labelsize=26)
    y1 = np.max(h[0])/len(VstdPC)
    if y1 % 0.1 >= 0.05:
        yf = 0.1*(int(y1*10) + 1)
    else:
        yf = 0.1*int(y1*10)

    yarr1 = np.arange(0, yf+0.01, 0.1)
    yarr0 = yarr1*len(VstdPC)

    plt.gca()
    plt.yticks(yarr0, np.round(yarr1, 1), fontsize=26)

    ax.yaxis.set_tick_params(labelsize=26)
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.text(0.15*x0 + 0.85*x1, 0.05*y0 + 0.95*y1, '(b)', fontsize=26)
    fig.savefig(path+'\\Figures\\Fig05b.png')
   
    
    
def fig06(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
        
    ISIfigures(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, path+'\\Figures\\Fig06.png')
    
def fig07(C_lags, autoC_mean, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
        
    autoCfigures(C_lags, autoC_mean, path+'\\Figures\\Fig07.png')
    
def fig08(frequency, power, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    SPDfigures(frequency, power, path+'\\Figures\\Fig08.png')



def fig09(cortex, xlim, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    Tools.raster_plot(cortex, xlim=xlim, savefig=path+'\\Figures\\Fig09.png')


def fig10(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    ISIfigures(ISImean, ISICV, C_lags, crossC_mean, Correlation_sigma, path+'\\Figures\\Fig10.png')

def fig11(C_lags, autoC_mean, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    autoCfigures(C_lags, autoC_mean, path+'\\Figures\\fig11.png')

def fig12(frequency, power, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    SPDfigures(frequency, power, path+'\\Figures\\fig12.png')
    
    
def fig13(cortex, pulse0, pulse1, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    
    p0_t0, p0_t1 = pulse0
    p1_t0, p1_t1 = pulse1
    Tools.raster_plot(cortex, xlim=(p0_t0-25, p0_t0+60), show=False)
    plt.vlines(p0_t0, 0, cortex.neurons.N+15, color='black', linestyle='dotted', linewidth=2)
    plt.vlines(p0_t0, min(cortex.neuron_idcs(('PC_L23', 0))), max(cortex.neuron_idcs(('PC_L23', 0))), color='purple', linestyle='--', linewidth=3)
    plt.text(p0_t0-38, 900, '(a)', fontsize=26)
    plt.savefig(path + '\\Figures\\Fig13a.png')
    
    Tools.raster_plot(cortex, xlim=(p1_t0-25, p1_t0+60), show=False)
    plt.vlines(p1_t0, 0, cortex.neurons.N+15, color='black', linestyle='dotted', linewidth=2)
    plt.vlines(p1_t0, min(cortex.neuron_idcs(('PC_L23', 0))), max(cortex.neuron_idcs(('PC_L23', 0))), color='purple', linestyle='--', linewidth=3)
    
    plt.text(p1_t0-38, 900, '(b)', fontsize=26)
    plt.savefig(path+'\\Figures\\Fig13b.png')
    
def fig14(cortex, pulse, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    p_t0, p_t1 = pulse
    
    Tools.raster_plot(cortex, xlim=(p_t0-25, p_t1+60), show=False)
    plt.vlines(p_t0, 0, cortex.neurons.N+15, color='black', linestyle='dotted', linewidth=2)
    plt.vlines(p_t0, min(cortex.neuron_idcs(('PC_L23', 0))), max(cortex.neuron_idcs(('PC_L23', 0))), color='purple', linestyle='--', linewidth=3)
    plt.text(p_t0-38, 900, '(a)', fontsize=26)
    plt.savefig(path + '\\Figures\\Fig14.png')
    
    
def fig15(Membr_std_list, PC_L23, PC_L5, path):
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    
    PC_L23_mean = []
    PC_L23_std = []
    PC_L5_mean = []
    PC_L5_std = []

    for scale in range(len(Membr_std_list)):
        PC_L23_mean.append(np.mean(PC_L23[scale]))
        PC_L23_std.append(np.std(PC_L23[scale]))
        PC_L5_mean.append(np.mean(PC_L5[scale]))
        PC_L5_std.append(np.std(PC_L5[scale]))
       
    
    variability =     Membr_std_list*100 
    PC_L23_mean = 100*np.asarray(PC_L23_mean)/PC_L23_mean[-1]
    PC_L5_mean = 100*np.asarray(PC_L5_mean)/PC_L5_mean[-1]

    PC_L23_SEM = 100*np.asarray(PC_L23_std)/np.sqrt(len(PC_L23[0]))/PC_L23_mean[-1]
    PC_L5_SEM = 100*np.asarray(PC_L5_std)/np.sqrt(len(PC_L5[0]))/PC_L5_mean[-1]

    
    plt.figure(figsize=(18, 12))
    plt.ylabel('relative spiking activity', fontsize=26)
    plt.xlabel('fraction of original membrane parameter SD (%)', fontsize=26)
    plt.tick_params(labelsize=26)
    plt.plot(variability, PC_L23_mean, label='PC L23', color='blue')
    plt.errorbar(variability, PC_L23_mean, yerr = PC_L23_SEM, fmt = 'o',color = 'blue', 
                ecolor = 'blue', elinewidth = 2, capsize=3)
    plt.plot(variability, PC_L5_mean, label='PC 5', color='orange')
    plt.errorbar(variability, PC_L5_mean, yerr = PC_L5_SEM, fmt = 'o',color = 'orange', 
                ecolor = 'orange', elinewidth = 2, capsize=3)
    plt.xlim(-20, 120)
    plt.xticks([0, 50, 100])
    plt.legend(prop={'size': 26})
    plt.savefig(path+'\\Figures\\Fig15.png')
    

def poissonfigures(cortex, pulse0, pulse1, file0, file1):
    PC_L23 = cortex.neuron_idcs(('PC_L23',0))
    PC_L5 = cortex.neuron_idcs(('PC_L5',0))
    
    p0_t0, p0_t1 = pulse0
    p1_t0, p1_t1 = pulse1
    
    Tools.raster_plot(cortex, xlim=(p0_t0-100, p0_t0+400), show=False)
    
    plt.vlines(p0_t0, 0, cortex.neurons.N+15, color='black', linestyle='dotted', linewidth=2)
    plt.vlines(p0_t1, 0, cortex.neurons.N+15, color='black', linestyle='dotted', linewidth=2)
    
    plt.vlines(p0_t0, min(PC_L23), max(PC_L23), color='purple', linestyle='--',  linewidth=3)
    plt.vlines(p0_t1, min(PC_L23), max(PC_L23), color='green', linestyle='--',  linewidth=3)
    
    plt.savefig(file0)
    
    
    Tools.raster_plot(cortex, xlim=(p1_t0-100, p1_t1+400), show=False)
    
    plt.vlines(p1_t0, 0, cortex.neurons.N+15, color='black', linestyle='dotted', linewidth=2)
    plt.vlines(p1_t1, 0, cortex.neurons.N+15, color='black', linestyle='dotted', linewidth=2)
    
    plt.vlines(p1_t0, min(PC_L5), max(PC_L5), color='purple', linestyle='--' , linewidth=3)
    plt.vlines(p1_t1, min(PC_L5), max(PC_L5), color='green', linestyle='--' , linewidth=3)
    
    plt.savefig(file1)


def fig16(cortex, pulse0, pulse1, path):    
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    poissonfigures(cortex, pulse0, pulse1, path+'\\Figures\\Fig16a.png', path+'\\Figures\\Fig16b.png')
    
def fig17(cortex, pulse0, pulse1, path):    
    if not os.path.isdir(path+'\\Figures'):
        os.mkdir(path+'\\Figures')
    poissonfigures(cortex, pulse0, pulse1, path+'\\Figures\\Fig17a.png', path+'\\Figures\\Fig17b.png')
