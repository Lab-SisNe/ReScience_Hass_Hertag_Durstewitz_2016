from matplotlib import pyplot as plt
import numpy as np, pandas as pd, brian2 as br2
from Auxiliar import time_report
from scipy.signal import periodogram
from scipy.stats import ttest_ind as ttest, mannwhitneyu as mwtest, chi2_contingency as chi2
from collections import namedtuple
from scipy.ndimage import gaussian_filter1d as gf1d


def raster_plot_simple(cortex, xlim=None, figsize=(18,10), s=3,
                       fontsize=24, labelsize=20, savefig=None, show=True):
    
    plt.figure(figsize=figsize)
    plt.scatter(cortex.spikemonitor.t_*1000, cortex.spikemonitor.i, s=s)
    
    plt.xlabel('time (ms)', fontsize=fontsize)
    plt.ylabel('neuron index', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    if xlim is not None:
        plt.xlim(xlim)
    
    if savefig is not None:
        plt.savefig(savefig)
     
    if show:
        plt.show()
        plt.close()
    
def raster_plot(cortex, xlim=None, figsize=(18,10), s=5, layerpadding=15,\
                fontsize=24, labelsize=20, lw=3, savefig=None, show=True):
    
    sp_t = cortex.spikemonitor.t_*1000
    sp_i = cortex.spikemonitor.i
    
    PC_L23idc = cortex.neuron_idcs(('PC_L23',0))
    PC_L23_isin = np.isin(sp_i, PC_L23idc)
    PC_L23_i = sp_i[PC_L23_isin]
    PC_L23_t = sp_t[PC_L23_isin]
    
    PC_L5idc = cortex.neuron_idcs(('PC_L5',0))
    PC_L5_isin = np.isin(sp_i, PC_L5idc)
    PC_L5_i = sp_i[PC_L5_isin]
    PC_L5_t = sp_t[PC_L5_isin]
    
    IN_L23idc = cortex.neuron_idcs(('IN_L23',0))
    IN_L23_isin = np.isin(sp_i, IN_L23idc)
    IN_L23_i = sp_i[IN_L23_isin]
    IN_L23_t = sp_t[IN_L23_isin]
    
    IN_L5idc = cortex.neuron_idcs(('IN_L5',0))
    IN_L5_isin = np.isin(sp_i, IN_L5idc)
    IN_L5_i = sp_i[IN_L5_isin]
    IN_L5_t = sp_t[IN_L5_isin]
    
    
    plt.figure(figsize=figsize)
    plt.scatter(PC_L23_t, PC_L23_i, s=s, color='blue')
    plt.scatter(IN_L23_t, IN_L23_i, s=s, color='red')
    plt.scatter(PC_L5_t, PC_L5_i+layerpadding, s=s, color='blue')
    plt.scatter(IN_L5_t, IN_L5_i+layerpadding, s=s, color='red')
 
    if xlim is None:
        xlim = cortex.transient, cortex.neurons.t/br2.ms
    
    plt.xlim(xlim)

    x0, x1 = plt.gca().get_xlim()
    y0, y1 = plt.gca().get_ylim()
    division_idc = cortex.network.group_distr[0][7][0] + layerpadding/2
    plt.hlines(division_idc, x0, x1, color='black', lw=lw)

    plt.xlabel('time (ms)', fontsize=fontsize)
    plt.ylabel('neuron index', fontsize=fontsize)
    
    plt.tick_params(labelsize=labelsize)

    yticks_space = cortex.network.basics.struct.Ncells//100*20
    yticks = np.arange(5*yticks_space+1)
    ytickspadding = yticks.copy()
    ytickspadding[yticks>division_idc] = ytickspadding[yticks>division_idc]+layerpadding

    if savefig is not None:
        plt.savefig(savefig)
    
    if show:
        plt.show()
        plt.close()

@time_report()
def get_V_stats(cortex, neuron_idcs=None, tlim=None, file=None):
    
    def extract_spikes(Varr, V_T, spike_jump=1):
        
        def get_postspike_positions(Varr, spike_jump=1): 
    
            Varr_diff = np.diff(np.abs(Varr))
            sppos = np.where(Varr_diff >= spike_jump)[0]
            
            return sppos+1
    
        def get_interspikes(Varr, post_spikes):
            interspikes = []
            i0 = 0
            for i1 in post_spikes:
                interspikes.append(Varr[i0:i1])
                i0=i1
            else:
                interspikes.append(Varr[i0:])
            return interspikes
    
        def extract_prespike(interspikes, V_T):
            for i in range(len(interspikes)-1):
                try:
                    V_T_crossing = np.max(np.where(interspikes[i]<V_T)[0])
                    interspikes[i] = (interspikes[i][:V_T_crossing+1])
                except ValueError:
                    interspikes[i] = np.asarray([])
    
        post_spikes = get_postspike_positions(Varr, spike_jump=1)
        interspikes = get_interspikes(Varr, post_spikes)
        extract_prespike(interspikes, V_T)
        
        return np.concatenate(interspikes)
    
    V = cortex.recorded.V.V/br2.mV

    if neuron_idcs is not None:
        V = V[neuron_idcs]
    else:
        neuron_idcs = np.arange(cortex.neurons.N)
    
    if tlim is not None:
        t0,t1=tlim
        t = cortex.recorded.V.t/br2.ms
        V = V[:, (t>=t0) & (t<t1)]
    
    Vmean = []
    Vstd = []
    V_minus_V_T_mean = []
    Vmean_extracted = []
    Vstd_extracted = []
    V_minus_V_T_mean_extracted = []
    
    
    for i in range(V.shape[0]):
        Varr = V[i]
        V_T = cortex.network.membr_params.loc[dict(param='V_T', cell_index=neuron_idcs[i])].values
        V_extracted = extract_spikes(Varr, V_T)
        
        Vmean.append(np.mean(Varr))
        Vstd.append(np.std(Varr))
        V_minus_V_T_mean.append(np.mean(Varr-V_T))
        Vmean_extracted.append(np.mean(V_extracted))   
        Vstd_extracted.append(np.std(V_extracted))
        V_minus_V_T_mean_extracted.append(np.mean(V_extracted - V_T))
            
    if file is not None:
        with open(file, 'w') as f:
            print('V correlations', file=f)
            print('{} cells'.format(len(neuron_idcs)), end='\n\n', file=f)
            print('Full V arrays', file=f)
            print('Vmean --> mean: {:.2f} mV, std: {:.2f} mV'.format(np.mean(Vmean), np.std(Vmean)), file=f)
            print('Vstd --> mean: {:.2f} mV, std: {:.2f} mV'.format(np.mean(Vstd), np.std(Vstd)), file=f)
            print('V minus V_T mean --> mean: {:.2f} mV, std: {:.2f} mV'.format(np.mean(V_minus_V_T_mean), np.std(V_minus_V_T_mean)), file=f, end='\n\n')
            print('V arrays after spike extraction', file=f)
            print('Vmean_extracted --> mean: {:.2f} mV, std: {:.2f} mV'.format(np.mean(Vmean_extracted), np.std(Vmean_extracted)), file=f)
            print('Vstd_extracted --> mean: {:.2f} mV, std: {:.2f} mV'.format(np.mean(Vstd_extracted), np.std(Vstd_extracted)), file=f)
            print('V minus V_T mean extracted--> mean: {:.2f} mV, std: {:.2f} mV'.format(np.mean(V_minus_V_T_mean_extracted), np.std(V_minus_V_T_mean_extracted)), file=f, end='\n\n')
    
    return Vmean,Vstd, V_minus_V_T_mean, Vmean_extracted, Vstd_extracted, V_minus_V_T_mean_extracted
  

@time_report('Correlations')
def get_correlations(cortex, tlim=None,  idcs=None, delta_t=2, lags=0, file=None, display=False, display_interval=10):
    
    def get_binned_spiking_trains(cortex, idcs, tlim, delta_t):
        if idcs is None:
            idcs = np.arange(cortex.neurons.N)
        elif isinstance(idcs, int):
            idcs = [idcs]
        
        if tlim is None:
            tlim = (cortex.transient, cortex.neurons.t/br2.ms)
        
        t0, t1 = tlim
        t1 = np.ceil(t1).astype(int)
        bin_arr = np.zeros((len(idcs), (t1-t0)//delta_t))
        
        spike_trains = cortex.spikemonitor.spike_trains()
        for i in range(len(idcs)):
            train = spike_trains[idcs[i]]/br2.ms
            train = train[(train>=t0)&(train<t1)]
            idc_arr =((train-t0)//delta_t).astype(int)
            bin_arr[i, idc_arr] = 1
            
        return bin_arr.astype(int)
        
    bin_arr = get_binned_spiking_trains(cortex, idcs, tlim, delta_t)
                       
    return ISIcorrelations(bin_arr, delta_t, lags, file, display, display_interval)

@time_report('original spike trains correlation')
def get_correlations_from_spike_trains(spike_trains, tlim, idcs=None, delta_t=2, lags=0, file=None, display=False, display_interval=10):
    
    def get_binned_spiking_trains(spike_trains, idcs, tlim, delta_t):
        if idcs is None:
            idcs = np.arange(len(spike_trains))
        elif isinstance(idcs, int):
            idcs = [idcs]
             
        t0,t1= tlim
        t1 = np.ceil(t1).astype(int)
        bin_arr = np.zeros((len(idcs), (t1-t0)//delta_t))
        
        for i in range(len(idcs)):
            train = spike_trains[idcs[i]]
            train = train[(train>=t0)&(train<t1)]
            idc_arr =((train-t0)//delta_t).astype(int)
            bin_arr[i, idc_arr] = 1
            
        return bin_arr.astype(int)
        
    bin_arr = get_binned_spiking_trains(spike_trains, idcs, tlim, delta_t)
                       
    return ISIcorrelations(bin_arr, delta_t, lags, file, display, display_interval)

def ISIcorrelations(bin_arr, delta_t, lags, file=None, display=False, display_interval=10):
     
    def group_correlations(bin_df, lag=0):
        
        def df_shifted(df, target, lag=0):
            cross = []
            auto = []
            for c in df.columns:       
                if c==target:
                    auto.append(df[target].corr(df[c].shift(periods=lag)))
                else:
                    cross.append(df[target].corr(df[c].shift(periods=lag)))

            return auto, cross
        
        auto = []
        cross = []
        for idc in range(bin_df.shape[1]):
            auto_idc, cross_idc = df_shifted(bin_df, idc, lag)
            auto.extend(auto_idc)
            cross.extend(cross_idc)
        return auto, cross

    if isinstance(lags, int):
        lags = [lags]
        

    bin_df = pd.DataFrame(bin_arr.transpose())
    auto_mean = []
    auto_std = []
    cross_mean = []
    cross_std = []
    curr = 0
    last = 0
    for l in lags:
        if display and curr==0 and last==0:
            print('Calculating correlations ...')
        auto, cross = group_correlations(bin_df, l)
        auto_mean.append(np.mean(auto))
        auto_std.append(np.std(auto))
        cross_mean.append(np.mean(cross))
        cross_std.append(np.std(cross))
        if display:
            curr += 100/len(lags)
            if curr>=display_interval:
                last += curr//display_interval
                curr=curr%display_interval
                print('{:.1f} % done'.format(last*display_interval+curr))
    
    lags = np.asarray(lags)*delta_t
    
    auto_mean = np.array(auto_mean)
    auto_std = np.asarray(auto_std)
    cross_mean = np.asarray(cross_mean)
    cross_std = np.asarray(cross_std)
    
    if file is not None:
        auto_zip = list(zip(lags, auto_mean, auto_std))
        cross_zip = list(zip(lags, cross_mean, cross_std))
        with open(file, 'w') as f:
            print('Auto-correlations:', file=f, end='\n\n')
            for lag in auto_zip:
                print('lag {} ms --> mean: {:.4f}, std: {:.4f}'.format(*lag), file=f)
            print(file=f)
            print('Auto-correlations:', file=f, end='\n\n')
            for lag in cross_zip:
                print('lag {} ms --> mean: {:.4f}, std: {:.4f}'.format(*lag), file=f)
                       
    return lags, np.array(auto_mean), np.array(auto_std), np.array(cross_mean), np.array(cross_std)

def show_correlations(bin_arr, lag, savefig=None, savedata=None):
    
    if savefig is not None:
        namesplit = savefig.split('.')
        autofig = namesplit.copy()
        autofig[-2] = autofig[-2]  + '_auto'
        autofig = '.'.join(autofig)
        crossfig = namesplit.copy()
        crossfig[-2] = crossfig[-2] + '_cross'
        crossfig='.'.join(crossfig)
         
    auto_mean, auto_std, cross_mean, cross_std = get_correlations(bin_arr, lag)

    lag_ms = lag*5
    
    plt.figure(figsize=(18,10))
    plt.plot(lag_ms, auto_mean, lw=3)
    plt.xlabel('lag (ms)', fontsize=24)
    plt.ylabel('autocorrelation', fontsize=24)
    plt.tick_params(labelsize=20)
    if savefig is not None:
        plt.savefig(autofig)
    plt.show()
    plt.close()


    plt.figure(figsize=(18,10))
    plt.plot(lag_ms, cross_mean, lw=3)
    plt.xlabel('lag (ms)', fontsize=24)
    plt.ylabel('crosscorrelation', fontsize=24)
    plt.tick_params(labelsize=20)
    if savefig is not None:
        plt.savefig(crossfig)
    plt.show()
    plt.close()
    
    corr_data = np.array([auto_mean, auto_std, cross_mean, cross_std]).transpose()
    corr_df = pd.DataFrame(corr_data, index=lag_ms, columns=['auto_mean', 'auto_std', 'cross_mean', 'cross_std'])
    corr_df.index.name = 'lag (ms)'
    
    if savedata is not None:
        with open(savedata, 'w') as f:
            print(corr_df, file=f)
            
    return corr_df
                  
def get_LFP(cortex, invert_Itot=True, population_agroupate=None):
    LFP = cortex.recorded.var('I_tot')/br2.pA
    t_LFP = cortex.recorded.t('I_tot')/br2.ms
    
    if population_agroupate is not None:
        if population_agroupate=='mean':
            LFP = np.mean(LFP, axis=0)
        elif population_agroupate=='sum':
            LFP = np.sum(LFP, axis=0)
        
    if invert_Itot:
        LFP = -LFP
        
    return t_LFP, LFP

def get_LFP_SPD(cortex, log=True, sigma=0):
    
    t, LFP = get_LFP(cortex, invert_Itot=False)
    fq = 1000/(t[1]-t[0])
    
    frequency, power = periodogram(LFP, fq)
    
    if log:
        power = np.log(power[frequency>0])
        frequency = np.log(frequency[frequency>0])
       
    
    if sigma>0:
        power = gf1d(power, sigma)
    
    return frequency, power
    
def get_LFP_SPD_from_Itotarray(Itotarray, fq, log=True, sigma=0):

    frequency, power = periodogram(Itotarray, fq)
    
    if log:
        power = np.log(power[frequency>0])
        frequency = np.log(frequency[frequency>0])
       
    
    if sigma>0:
        power = gf1d(power, sigma)
    
    return frequency, power
    
    

# def get_ISI_stats(cortex, rate=None, groupstripe=None, tlim=None, savetxt=None):
def get_ISI_stats(cortex, neuron_idcs=None, tlim=None, savetxt=None):
    
    spike_trains = [*cortex.spikemonitor.spike_trains().values()]
    
    if neuron_idcs is None:
        neuron_idcs = np.arange(spike_trains)
    spike_trains = [spike_trains[idc]/br2.ms for idc in neuron_idcs]
    
    if tlim is not None:
        t0, t1 = tlim
        spike_trains = [train[(train>=t0)&(train<t1)] for train in spike_trains]
     
    return ISIstats(spike_trains, savetxt)
   
 
def get_ISI_stats_from_spike_trains(spike_trains, neuron_idcs=None, tlim=None, savetxt=None):
       
    if neuron_idcs is None:
        neuron_idcs = np.arange(len(spike_trains))
    spike_trains = [spike_trains[idc] for idc in neuron_idcs]
    
    
    if tlim is not None:
        t0, t1 = tlim
        spike_trains = [train[(train>=t0)&(train<t1)] for train in spike_trains]
        
    return ISIstats(spike_trains, savetxt)


def ISIstats(spike_trains, savetxt):
    
    ISImean = []
    ISICV = []
    
    for train in spike_trains:
        ISI = np.diff(train)
        ISImean.append(np.mean(ISI))
        ISICV.append(np.std(ISI)/np.mean(ISI))
    
    if savetxt is not None:
        with open(savetxt,'w') as f:
            print('ISI stats\n', file=f)
            print('Cells:', len(ISImean), file=f)
            print('ISImean mean: {:.3f} ms'.format(np.mean(ISImean)), file=f)
            print('ISImean std: {:.3f}ms'.format(np.std(ISImean)), file=f)
            print('ISICV mean: {:.3f}'.format(np.mean(ISICV)), file=f)
            print('ISICV std: {:.3f}'.format(np.std(ISICV)), file=f)
        
    return ISImean, ISICV



def comp_membrparam_rategroup(cortex, rate, groupstripe_list, file=None):
    
    def save_membrparam_rategroup(group_dict, rate, file):
        
        with open(file, 'w') as f:
            for groupstripe in group_dict:
                print(groupstripe, file=f, end='\n\n')

                for param in group_dict[groupstripe]:
                    less = group_dict[groupstripe][param]['less']
                    geq = group_dict[groupstripe][param]['greater_equal']
                    mw = group_dict[groupstripe][param]['mwtest']
                    tt = group_dict[groupstripe][param]['ttest']
                    
                    print('{}:'.format(param), file=f)
                    print('rate <  {} Hz --> mean: {:.2f}, std: {:.2f} ({} cells)'.format(rate, less.mean, less.std, less.N), file=f)
                    print('rate >= {} Hz --> mean: {:.2f}, std: {:.2f} ({} cells)'.format(rate, geq.mean, geq.std, geq.N), file=f, end='\n\n')
                    
                    if mw.pvalue>=0.05:
                        mw_star = ''
                    elif mw.pvalue>=0.01:
                        mw_star = '*'
                    elif mw.pvalue >= 0.001:
                        mw_star = '**'
                    else:
                        mw_star = '***'
                    
                    if tt.pvalue>=0.05:
                        tt_star = ''
                    elif tt.pvalue>=0.01:
                        tt_star = '*'
                    elif tt.pvalue >= 0.001:
                        tt_star = '**'
                    else:
                        tt_star = '***'
                    
                    print('Mann-Whitney U: {:.2f}, p-value: {:.3f}'.format(mw.statistic, mw.pvalue)+mw_star, file=f)
                    print("Student's t: {:.2f}, p-value: {:.3f}".format(tt.statistic, tt.pvalue)+tt_star, file=f)
                    print('-'*40, file=f)
                print('='*40, file=f)
                print('='*40, file=f)
             
    
    Stats = namedtuple('Stats', ['mean', 'std', 'N'])
    
    if isinstance(groupstripe_list[0], (int, str)):
        groupstripe_list = [groupstripe_list]
         
    # less_group_dict = {}
    # geq_group_dict = {}
    # mwtest_group_dict = {}
    
    group_dict = {}
    
    for groupstripe in groupstripe_list:
        group, stripe = groupstripe
        gsname = '{}_stripe_{}'.format(group, stripe)
        
        cell_less = cortex.spiking_idcs((np.less, rate), groupstripe)
        cell_geq = cortex.spiking_idcs((np.greater_equal, rate), groupstripe)
        
        param_dict = {}
        
        for param in cortex.network.membr_params.coords['param'].values:
            less_params = cortex.network.membr_params.loc[dict(param=param, cell_index=cell_less)].values
            geq_params = cortex.network.membr_params.loc[dict(param=param, cell_index=cell_geq)].values

            param_dict[param] = {}
            param_dict[param]['less'] = Stats(np.mean(less_params), np.std(less_params), len(less_params))
            param_dict[param]['greater_equal'] = Stats(np.mean(geq_params), np.std(geq_params), len(geq_params))
            param_dict[param]['mwtest'] = mwtest(less_params, geq_params, alternative='two-sided')
            param_dict[param]['ttest'] = ttest(less_params, geq_params)
   
        rheo_less = cortex.rheobase(cell_less) 
        rheo_geq = cortex.rheobase(cell_geq)
        
        param_dict['Rheobase'] = {}
        param_dict['Rheobase']['less'] = Stats(np.mean(rheo_less), np.std(rheo_less), len(rheo_less))
        param_dict['Rheobase']['greater_equal'] = Stats(np.mean(rheo_geq), np.std(rheo_geq), len(rheo_geq))
        param_dict['Rheobase']['mwtest'] = mwtest(rheo_less, rheo_geq, alternative='two-sided')
        param_dict['Rheobase']['ttest'] = ttest(rheo_less, rheo_geq)
    
        group_dict[gsname] = param_dict
     
    if file is not None:
        save_membrparam_rategroup(group_dict, rate, file)
        
    return group_dict

def contingency(cortex, rate, target_groupstripe, source_groupstripe, channel=None, file=None):
    
    def save_contingency(group_dict, rate, file, channel):
        
        with open(file, 'w') as f:
            for grouptarget in group_dict:
                cont_tab, N_less, N_geq, pCon_less, pCon_geq, chi2_res, pvalue = group_dict[grouptarget]
                if channel is not None:
                    print('Channel:', channel, file=f, end='\n\n')
                print(grouptarget, file=f, end='\n\n')
                print('rate <  {} Hz --> pCon = {:.2f}% ({} cells)'.format(rate, pCon_less*100, N_less), file=f)
                print('rate >= {} Hz --> pCon = {:.2f}% ({} cells)'.format(rate, pCon_geq*100, N_geq), file=f, end='\n\n')
                
                print('Contingency table', file=f)
                print('[{}]'.format(cont_tab[0]), file=f)
                print('[{}]'.format(cont_tab[1]), file=f, end='\n\n')
                if pvalue >= 0.05:
                    chistar = ''
                elif pvalue >= 0.01:
                    chistar = '*'
                elif pvalue >= 0.001:
                    chistar = '**'
                else:
                    chistar = '***'
                print('chi2 = {:.2f}, p-value: {:.3f}'.format(chi2_res, pvalue)+chistar, file=f)
                print('-'*40, file=f, end='\n\n')
                
    if isinstance(target_groupstripe[0], (int, str)):
        target_groupstripe = [target_groupstripe]
         
    if isinstance(source_groupstripe[0], (int, str)):
        source_groupstripe = [source_groupstripe]
    
    group_dict = {}
    
    for source in source_groupstripe:
        source_group, source_stripe = source
        source_cell = cortex.neuron_idcs(source_groupstripe)
    
        for target in target_groupstripe:
            target_group, target_stripe = target  
            target_cell_less = cortex.spiking_idcs((np.less, rate), target)
            target_cell_geq = cortex.spiking_idcs((np.greater_equal, rate), target)
           
            gsname = 'to_{}_stripe_{}_from_{}_stripe_{}'.format(target_group, target_stripe, source_group, source_stripe)
          
        
            
            
            
            syn_less = cortex.syn_idcs_from_neurons(target_cell_less, source_cell, channel)
            syn_geq = cortex.syn_idcs_from_neurons(target_cell_geq, source_cell, channel)
            
            Ncon_less = len(syn_less)
            Ntot_less = len(target_cell_less) * len(source_cell)
            Ndis_less = Ntot_less - Ncon_less


            
            Ncon_geq = len(syn_geq)
            Ntot_geq = len(target_cell_geq) * len(source_cell)
            Ndis_geq = Ntot_geq - Ncon_geq
            
            contingency_tab = [[Ncon_less, Ndis_less], [Ncon_geq, Ndis_geq]]
            chi2_result, pvalue = chi2(contingency_tab)[:2]
            
            group_dict[gsname] = contingency_tab, len(target_cell_less), len(target_cell_geq),\
                                Ncon_less/Ntot_less, Ncon_geq/Ntot_geq, chi2_result, pvalue
     
    if file is not None:
        save_contingency(group_dict, rate, file, channel)
        
    return group_dict
    
def comp_synparam_rategroup(cortex, rate, target_groupstripe, source_groupstripe, channel=None, file=None):
    
    def save_synparam_rategroup(group_dict, rate, file, channel):
        with open(file, 'w') as f:
        
            for grouptarget in group_dict:
                if channel is not None:
                    print('Channel:', channel, file=f, end='\n\n')
                
                print(grouptarget, file=f, end='\n\n')
                # print(grouptarget)
                # print(group_dict[grouptarget])
                for param in group_dict[grouptarget]:
                    less = group_dict[grouptarget][param]['less']
                    geq = group_dict[grouptarget][param]['greater_equal']
                    mw = group_dict[grouptarget][param]['mwtest']
                    tt = group_dict[grouptarget][param]['ttest']
                    
                    print('{}:'.format(param), file=f)
                    print('rate <  {} Hz --> mean: {:.2f}, std: {:.2f} ({} synapses)'.format(rate, less.mean, less.std, less.N), file=f)
                    print('rate >= {} Hz --> mean: {:.2f}, std: {:.2f} ({} synapses)'.format(rate, geq.mean, geq.std, geq.N), file=f, end='\n\n')
                    
                    if mw.pvalue>=0.05:
                        mw_star = ''
                    elif mw.pvalue>=0.01:
                        mw_star = '*'
                    elif mw.pvalue >= 0.001:
                        mw_star = '**'
                    else:
                        mw_star = '***'
                    
                    if tt.pvalue>=0.05:
                        tt_star = ''
                    elif tt.pvalue>=0.01:
                        tt_star = '*'
                    elif tt.pvalue >= 0.001:
                        tt_star = '**'
                    else:
                        tt_star = '***'
                    
                    print('Mann-Whitney U: {:.2f}, p-value: {:.3f}'.format(mw.statistic, mw.pvalue)+mw_star, file=f)
                    print("Student's t: {:.2f}, p-value: {:.3f}".format(tt.statistic, tt.pvalue)+tt_star, file=f)
                    print('-'*40, file=f)
                print('='*40, file=f)
                print('='*40, file=f)
             
            
       
    if isinstance(target_groupstripe[0], (int, str)):
        target_groupstripe = [target_groupstripe]
         
    if isinstance(source_groupstripe[0], (int, str)):
        source_groupstripe = [source_groupstripe]
    
    Stats = namedtuple('Stats', ['mean', 'std', 'N'])
    group_dict = {}
    
    for source in source_groupstripe:
        source_group, source_stripe = source
        source_cell = cortex.neuron_idcs(source_groupstripe)
    
        for target in target_groupstripe:
            target_group, target_stripe = target  
            target_cell_less = cortex.spiking_idcs((np.less, rate), target)
            target_cell_geq = cortex.spiking_idcs((np.greater_equal, rate), target)
           
            gsname = 'to_{}_stripe_{}_from_{}_stripe_{}'.format(target_group, target_stripe, source_group, source_stripe)
            
            syn_less = cortex.syn_idcs_from_neurons(target_cell_less, source_cell, channel)
            syn_geq = cortex.syn_idcs_from_neurons(target_cell_geq, source_cell, channel)
            # print(syn_less)
            # input()
            # return syn_less
            param_dict={}
            for k in cortex.network.syn_params:
                for param in cortex.network.syn_params[k].coords['param'].values:
                    less_params = cortex.network.syn_params[k].loc[dict(param=param, syn_index=syn_less)].values
                   
                    geq_params = cortex.network.syn_params[k].loc[dict(param=param, syn_index=syn_geq)].values

                    param_dict[param] = {}
                    param_dict[param]['less'] = Stats(np.mean(less_params), np.std(less_params), len(less_params))
                    param_dict[param]['greater_equal'] = Stats(np.mean(geq_params), np.std(geq_params), len(geq_params))
                    param_dict[param]['mwtest'] = mwtest(less_params, geq_params, alternative='two-sided')
                    param_dict[param]['ttest'] = ttest(less_params, geq_params)
                      
            
            group_dict[gsname] = param_dict

            

    if file is not None:
        save_synparam_rategroup(group_dict, rate, file, channel)
        
    return group_dict
    
def get_spiking(cortex, rate, groupstripe, file=None):
    Ntotal = cortex.network.basics.struct.Ncells_total
    Nspiking = len(cortex.spiking_idcs((np.greater_equal, rate), groupstripe))
    Nnotspiking = len(cortex.spiking_idcs((np.less, rate), groupstripe))
    
    Pspiking = Nspiking/Ntotal
    Pnotspiking = Nnotspiking/Ntotal
    
    if file is not None:
        with open(file, 'w') as f:
            print('Group_stripe:', groupstripe, file=f)
            print('Spiking rate: >= {} Hz'.format(rate), file=f)
            print('Total N: {} cells'.format(Ntotal), file=f, end='\n\n')
            print('Spiking: {:.2f} % ({} cells)'.format(Pspiking*100, Nspiking), file=f)
            print('Not spiking: {:.2f} % ({} cells)'.format(Pnotspiking*100, Nnotspiking), file=f)
            
    
    return Pspiking, Pnotspiking


def get_membr_params(cortex, groupstripe_list, alias_list=None, file=None):
    
    
    if isinstance(groupstripe_list[0], (str, int)):
        groupstripe_list = [groupstripe_list]
    
    group_dict = {}
    for groupstripe in groupstripe_list:
        gsname = 'group_{}_stripe_{}'.format(*groupstripe)
        neuron_idcs = cortex.neuron_idcs(groupstripe)
        
        param_dict = {}
        for param in cortex.network.membr_params.coords['param'].values:
            
            
            param_values = cortex.network.membr_params.loc[dict(param=param, cell_index=neuron_idcs)].values
            param_dict[param] = np.mean(param_values), np.std(param_values)
        
        group_dict[gsname] = param_dict
     
    if alias_list is not None:
        alias_dict={}
        for gs in range(len(group_dict)):
            alias_dict[list(group_dict.keys())[gs]] = alias_list[gs] 
     
    if file is not None:
        with open(file, 'w') as f:
            for gsname in group_dict:
                if alias_list is not None:
                    print(alias_dict[gsname], '({})'.format(gsname), file=f, end='\n\n')
                else:
                    print(gsname, end='\n\n', file=f)
                for param in group_dict[gsname]:
                    print(param, '--> mean: {:.2f}, std: {:.2f}'.format(*group_dict[gsname][param]), file=f)
               
                print('-'*30, file=f, end='\n\n')
            print('='*40, file=f, end='\n\n')
    
    return group_dict

if __name__ == '__main__':
    
    # aliases = ['PC', 'fast-spiking cells','bitufted cells', 'basket cells', 'Matinotti cells']
    # get_membr_params(cortex,[('PC',0), ('IN_L_both',0), ('IN_CL_both',0), ('IN_CC',0), ('IN_F',0)], alias_list = aliases, file='compWWW.txt')
    pass