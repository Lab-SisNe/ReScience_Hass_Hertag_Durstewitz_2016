import numpy as np
import xarray as xr
from dataclasses import dataclass, field
from numpy.matlib import repmat
from collections import namedtuple
import brian2 as br2
from Auxiliar import time_report


GROUP_KINDS = {
    'PC_L23': 'exc', 'IN_L_L23': 'inh', 'IN_L_d_L23': 'inh', 'IN_CL_L23': 'inh', 'IN_CL_AC_L23': 'inh', 'IN_CC_L23': 'inh', 'IN_F_L23': 'inh',
    'PC_L5': 'exc', 'IN_L_L5': 'inh', 'IN_L_d_L5': 'inh', 'IN_CL_L5': 'inh', 'IN_CL_AC_L5': 'inh', 'IN_CC_L5': 'inh', 'IN_F_L5': 'inh',
    }

GROUP_NAMES = list(GROUP_KINDS.keys())

MEMB_PARAM_UNITS = dict(C=dict(unit='pF', value='C'), g_L=dict(unit='nS', value='g_L'),
                        E_L=dict(unit='mV', value='E_L'), delta_T=dict(unit='mV', value='delta_T'),
                        V_up=dict(unit='mV', value='V_up'), tau_w=dict(unit='ms', value='tau_w'),
                        b=dict(unit='pA', value='b'), V_r=dict(unit='mV', value='V_r'),
                        V_T=dict(unit='mV', value='V_T'),
                        )

UNITS_DICT = {'pF':br2.pF, 'nS': br2.nS, 'mV': br2.mV,
              'ms': br2.ms, 'pA': br2.pA, 1:1,
              }


UNIT_MAIN_DICT = {'pF': 'farad', 'nS': 'siemens', 'mV': 'volt', 'ms': 'second',
                  'pA': 'amp', 1: 1}

MEMB_PARAM_NAMES = list(MEMB_PARAM_UNITS.keys())
MembraneTuple = namedtuple('MembraneTuple', MEMB_PARAM_NAMES
                                              )



CHANNEL_KINDS={'AMPA': 'exc', 'GABA':'inh', 'NMDA':'exc'}
CHANNEL_NAMES = [name for name in list(CHANNEL_KINDS.keys())]


CHANNEL_PARAMS= ['tau_on', 'tau_off', 'E', 'Mg_fac',
                 'Mg_slope', 'Mg_half']
    
CHANNEL_PARAMS_UNITS = dict(tau_on=('ms','second'), tau_off=('ms', 'second'),
                            E=('mV', 'volt'), Mg_fac=('1','1'),
                            Mg_slope=('1','1'), Mg_half=('1','1'))

PCELLS_PER_GROUP = np.asarray([47, 1.55, 1.55, 1.3, 1.3, 2.6, 2.1,
                                  38, 0.25, 0.25, 0.25, 0.25, 1.8, 1.8 ])



GROUPS_PERTAINING={
    'PC_L23':['PC_L23', 'PC', 'L23', 'ALL'], 
    'IN_L_L23':['IN_L_L23', 'IN_L23', 'IN_L', 'IN', 'L23', 'ALL'], 
    'IN_L_d_L23':['IN_L_d_L23','IN_L23', 'IN_L_d', 'IN','L23', 'ALL'], 
    'IN_CL_L23':['IN_CL_L23','IN_L23', 'IN_CL', 'IN', 'L23', 'ALL'],
    'IN_CL_AC_L23':['IN_CL_AC_L23','IN_L23', 'IN_CL_AC', 'IN', 'L23', 'ALL'], 
    'IN_CC_L23':['IN_CC_L23','IN_L23', 'IN_CC', 'IN', 'L23', 'ALL'], 
    'IN_F_L23':['IN_F_L23','IN_L23', 'IN_F', 'IN', 'L23', 'ALL'],
    'PC_L5':['PC_L5', 'PC', 'L5', 'ALL'], 
    'IN_L_L5':['IN_L_L5', 'IN_L5', 'IN_L', 'IN', 'L5', 'ALL'], 
    'IN_L_d_L5':['IN_L_d_L5','IN_L5', 'IN_L_d', 'IN', 'L5', 'ALL'], 
    'IN_CL_L5':['IN_CL_L5','IN_L5', 'IN_CL', 'IN', 'L5', 'ALL'],
    'IN_CL_AC_L5':['IN_CL_AC_L5','IN_L5', 'IN_CL_AC', 'IN', 'L5', 'ALL'], 
    'IN_CC_L5':['IN_CC_L5','IN_L5', 'IN_CC', 'IN', 'L5', 'ALL'], 
    'IN_F_L5':['IN_F_L5','IN_L5', 'IN_F', 'IN', 'L5', 'ALL'],
    }


GROUPS_SETS={}
for group in GROUPS_PERTAINING:
    for gr in GROUPS_PERTAINING[group]:
        if gr in GROUPS_SETS:
            GROUPS_SETS[gr].append(group)
        else:
            GROUPS_SETS[gr]=[group,]


GROUPS_SETS = {
    'ALL': ['PC_L23','IN_L_L23','IN_L_d_L23','IN_CL_L23','IN_CL_AC_L23','IN_CC_L23','IN_F_L23',
            'PC_L5','IN_L_L5','IN_L_d_L5','IN_CL_L5','IN_CL_AC_L5','IN_CC_L5','IN_F_L5'],
    
    'L23': ['PC_L23','IN_L_L23','IN_L_d_L23','IN_CL_L23','IN_CL_AC_L23','IN_CC_L23','IN_F_L23'],
    'L5': ['PC_L5', 'IN_L_L5','IN_L_d_L5','IN_CL_L5','IN_CL_AC_L5','IN_CC_L5','IN_F_L5'],
    
    'PC': ['PC_L23', 'PC_L5'], 
    'PC_L23': ['PC_L23'],
    'PC_L5': ['PC_L5'],
    
    'IN': ['IN_L_L23', 'IN_L_d_L23','IN_CL_L23','IN_CL_AC_L23','IN_CC_L23','IN_F_L23',
           'IN_L_L5','IN_L_d_L5','IN_CL_L5','IN_CL_AC_L5','IN_CC_L5','IN_F_L5'],
    'IN_L23': ['IN_L_L23', 'IN_L_d_L23','IN_CL_L23','IN_CL_AC_L23','IN_CC_L23','IN_F_L23'],
    'IN_L5': ['IN_L_L5', 'IN_L_d_L5','IN_CL_L5','IN_CL_AC_L5','IN_CC_L5','IN_F_L5'],
    
    'IN_L_both': ['IN_L_L23', 'IN_L_d_L23', 'IN_L_L5', 'IN_L_d_L5'],
    'IN_CL_both': ['IN_CL_L23', 'IN_CL_AC_L23', 'IN_CL_L5', 'IN_CL_AC_L5'],
    'IN_L_both_L23': ['IN_L_L23', 'IN_L_d_L23'],
    'IN_L_both_L5': ['IN_L_L5', 'IN_L_d_L5'],
    'IN_CL_both_L23': ['IN_CL_L23', 'IN_CL_AC_L23'],
    'IN_CL_both_L5': ['IN_CL_L5', 'IN_CL_AC_L5'],
    
    'IN_L': ['IN_L_L23', 'IN_L_L5'],
    'IN_L_d': ['IN_L_d_L23', 'IN_L_d_L5'],
    'IN_CL': ['IN_CL_L23', 'IN_CL_L5'],
    'IN_CL_AC': ['IN_CL_AC_L23', 'IN_CL_AC_L5'],
    'IN_CC': ['IN_CC_L23', 'IN_CC_L5'],
    'IN_F': ['IN_F_L23', 'IN_F_L5'],
     
    'IN_L_L23': ['IN_L_L23'],
    'IN_L_d_L23': ['IN_L_d_L23'],
    'IN_CL_L23': ['IN_CL_L23'],
    'IN_CL_AC_L23': ['IN_CL_AC_L23'],
    'IN_CC_L23': ['IN_CC_L23'],
    'IN_F_L23': ['IN_F_L23'],
    
    'IN_L_L5': ['IN_L_L5'],
    'IN_L_d_L5': ['IN_L_d_L5'],
    'IN_CL_L5': ['IN_CL_L5'],
    'IN_CL_AC_L5': ['IN_CL_AC_L5'],
    'IN_CC_L5': ['IN_CC_L5'],
    'IN_F_L5': ['IN_F_L5']}

CELL_MEAN = np.zeros((len(MEMB_PARAM_NAMES), len(GROUP_NAMES)))
CELL_MEAN = xr.DataArray(CELL_MEAN, coords=[MEMB_PARAM_NAMES, GROUP_NAMES], dims=['param', 'group'], name='Memb_param mean')

CELL_MEAN.loc[dict(param='C')] = [ 3.0751, 1.6902, 1.6902, 3.0014, 3.0014, 3.0751, 3.3869,
                                      2.2513, 1.6902, 1.6902, 3.0014, 3.0014, 2.2513, 3.3869, ]
CELL_MEAN.loc[dict(param='g_L')] = [ 1.9661, 1.0353, 1.0353, 1.4581, 1.4581, 1.9661, 1.0106, 
                                        1.0196, 1.0353, 1.0353, 1.4581, 1.4581, 1.0196, 1.0106, ]
CELL_MEAN.loc[dict(param='E_L')] = [ 3.5945, 2.9528, 2.9528, 3.0991, 3.0991, 3.5945, 3.8065,
                                        3.4415, 2.9528, 2.9528, 3.0991, 3.0991, 3.4415, 3.8065, ]
CELL_MEAN.loc[dict(param='delta_T')] = [ 1.0309, 3.2163, 3.2163, 3.1517, 3.1517, 1.0309, 3.0269, 
                                            1.5178, 3.2163, 3.2163, 3.1517, 3.1517, 1.5178, 3.0269, ]
CELL_MEAN.loc[dict(param='V_up')] = [ 3.1428, 2.8230, 2.8230, 2.9335, 2.9335, 3.1428, 2.3911, 
                                         1.0702, 2.8230, 2.8230, 2.9335, 2.9335, 1.0702, 2.3911, ]
CELL_MEAN.loc[dict(param='tau_w')] = [ 4.4809, 1.0542, 1.0542, 1.0730, 1.0730, 4.4809, 4.1986, 
                                          4.5650, 1.0542, 1.0542, 1.0730, 1.0730, 4.5650, 4.1986, ]
CELL_MEAN.loc[dict(param='b')] = [ 1.0189, 2.5959, 2.5959, 0.6931, 0.6931, 1.0189, 0.8080, 
                                      1.1154, 2.5959, 2.5959, 0.6931, 0.6931, 1.1154, 0.8080, ]
CELL_MEAN.loc[dict(param='V_r')] = [ 5.0719, 4.1321, 4.1321, 1.9059, 1.9059, 5.0719, 3.0051, 
                                        4.3414, 4.1321, 4.1321, 1.9059, 1.9059, 4.3414, 3.0051, ]
CELL_MEAN.loc[dict(param='V_T')] = [ 2.9010, 3.6925, 3.6925, 2.9462, 2.9462, 2.9010, 3.0701,
                                        3.3302, 3.6925, 3.6925, 2.9462, 2.9462, 3.3302, 3.0701, ]

CELL_COVARIANCE = np.zeros((len(GROUP_NAMES), len(MEMB_PARAM_NAMES), len(MEMB_PARAM_NAMES)))
CELL_COVARIANCE=xr.DataArray(CELL_COVARIANCE, coords=[GROUP_NAMES, MEMB_PARAM_NAMES,MEMB_PARAM_NAMES], dims=['group', 'memb_par0', 'memb_par1'], name='Memb_param covariance')
                            
CELL_COVARIANCE.loc[dict(group='PC_L23')] = np.matrix([[1.0000, 0.1580, -0.5835, 0.4011, -0.0561, 0.0718, -0.2038, 0.2615, -0.2365,],
                                [0.1580, 1.0000, 0.0141, -0.1272, -0.4327, 0.1778, -0.0902, -0.0329, -0.3778,],
                                [-0.5835, 0.0141, 1.0000, -0.6295, -0.2949, -0.2008, 0.3164, -0.2615, -0.0536,],
                                [0.4011, -0.1272, -0.6295, 1.0000, 0.6960, -0.2587, -0.0988, 0.6113, 0.5636,],
                                [-0.0561, -0.4327, -0.2949, 0.6960, 1.0000, -0.3370, 0.2042, 0.3959, 0.8581,],
                                [0.0718, 0.1778, -0.2008, -0.2587, -0.3370, 1.0000, -0.0634, -0.5202, -0.3829,],
                                [-0.2038, -0.0902, 0.3164, -0.0988, 0.2042, -0.0634, 1.0000, 0.0559, 0.3322,],
                                [0.2615, -0.0329, -0.2615, 0.6113, 0.3959, -0.5202, 0.0559, 1.0000, 0.3210,],
                                [-0.2365, -0.3778, -0.0536, 0.5636, 0.8581, -0.3829, 0.3322, 0.3210, 1.0000,], ])
CELL_COVARIANCE.loc[dict(group='IN_L_L23')] = np.matrix([[1.0000, -0.2894, 0.0381, 0.0664, -0.2418, 0.2253, 0.2822, -0.2919, 0.0581,],
                                                   [-0.2894, 1.0000, -0.2259, 0.4265, 0.1859, -0.6307, -0.0140, 0.4944, 0.2495,],
                                                   [0.0381, -0.2259, 1.0000, -0.2855, 0.0724, 0.1199, -0.1487, -0.3773, 0.1881,],
                                                   [0.0664, 0.4265, -0.2855, 1.0000, 0.2208, -0.3752, 0.0660, 0.3415, 0.7289,],
                                                   [-0.2418, 0.1859, 0.0724, 0.2208, 1.0000, 0.1412, -0.2931, 0.1993, 0.4609,],
                                                   [0.2253, -0.6307, 0.1199, -0.3752, 0.1412, 1.0000, -0.2855, -0.2046, -0.1974,],
                                                   [0.2822, -0.0140, -0.1487, 0.0660, -0.2931, -0.2855, 1.0000, -0.1172, -0.0851,],
                                                   [-0.2919, 0.4944, -0.3773, 0.3415, 0.1993, -0.2046, -0.1172, 1.0000, 0.0530,],
                                                   [0.0581, 0.2495, 0.1881, 0.7289, 0.4609, -0.1974, -0.0851, 0.0530, 1.0000], ])
CELL_COVARIANCE.loc[dict(group='IN_L_d_L23')] = np.matrix([[1.0000, -0.2894, 0.0381, 0.0664, -0.2418, 0.2253, 0.2822, -0.2919, 0.0581,],
                               [-0.2894, 1.0000, -0.2259, 0.4265, 0.1859, -0.6307, -0.0140, 0.4944, 0.2495,],
                               [0.0381, -0.2259, 1.0000, -0.2855, 0.0724, 0.1199, -0.1487, -0.3773, 0.1881,],
                               [0.0664, 0.4265, -0.2855, 1.0000, 0.2208, -0.3752, 0.0660, 0.3415, 0.7289,],
                               [-0.2418, 0.1859, 0.0724, 0.2208, 1.0000, 0.1412, -0.2931, 0.1993, 0.4609,],
                               [0.2253, -0.6307, 0.1199, -0.3752, 0.1412, 1.0000, -0.2855, -0.2046, -0.1974,],
                               [0.2822, -0.0140, -0.1487, 0.0660, -0.2931, -0.2855, 1.0000, -0.1172, -0.0851,],
                               [-0.2919, 0.4944, -0.3773, 0.3415, 0.1993, -0.2046, -0.1172, 1.0000, 0.0530,],
                               [0.0581, 0.2495, 0.1881, 0.7289, 0.4609, -0.1974, -0.0851, 0.0530, 1.0000], ])
CELL_COVARIANCE.loc[dict(group='IN_CL_L23')] = np.matrix([[1.0000, -0.2394, -0.6001, 0.3114, -0.2367, 0.5856, 0.2077, 0.0171, -0.4079,],
                              [-0.2394, 1.0000, -0.1764, 0.4675, 0.1810, -0.4942, -0.4389, 0.6950, 0.0811,],
                              [-0.6001, -0.1764, 1.0000, -0.6002, 0.2170, -0.0922, 0.2129, -0.3566, 0.4204,],
                              [0.3114, 0.4675, -0.6002, 1.0000, 0.2597, -0.1039, -0.5507, 0.7230, 0.0775,],
                              [-0.2367, 0.1810, 0.2170, 0.2597, 1.0000, 0.2159, -0.7123, 0.0193, 0.8494,],
                              [0.5856, -0.4942, -0.0922, -0.1039, 0.2159, 1.0000, 0.0587, -0.4724, 0.0957,],
                              [0.2077, -0.4389, 0.2129, -0.5507, -0.7123, 0.0587, 1.0000, -0.3395, -0.5780,],
                              [0.0171, 0.6950, -0.3566, 0.7230, 0.0193, -0.4724, -0.3395, 1.0000, -0.1084,],
                              [-0.4079, 0.0811, 0.4204, 0.0775, 0.8494, 0.0957, -0.5780, -0.1084, 1.0000], ])
CELL_COVARIANCE.loc[dict(group='IN_CL_AC_L23')] = np.matrix([[1.0000, -0.2394, -0.6001, 0.3114, -0.2367, 0.5856, 0.2077, 0.0171, -0.4079,],
                                 [-0.2394, 1.0000, -0.1764, 0.4675, 0.1810, -0.4942, -0.4389, 0.6950, 0.0811,],
                                 [-0.6001, -0.1764, 1.0000, -0.6002, 0.2170, -0.0922, 0.2129, -0.3566, 0.4204,],
                                 [0.3114, 0.4675, -0.6002, 1.0000, 0.2597, -0.1039, -0.5507, 0.7230, 0.0775,],
                                 [-0.2367, 0.1810, 0.2170, 0.2597, 1.0000, 0.2159, -0.7123, 0.0193, 0.8494,],
                                [0.5856, -0.4942, -0.0922, -0.1039, 0.2159, 1.0000, 0.0587, -0.4724, 0.0957,],
                                 [0.2077, -0.4389, 0.2129, -0.5507, -0.7123, 0.0587, 1.0000, -0.3395, -0.5780,],
                                 [0.0171, 0.6950, -0.3566, 0.7230, 0.0193, -0.4724, -0.3395, 1.0000, -0.1084,],
                                 [-0.4079, 0.0811, 0.4204, 0.0775, 0.8494, 0.0957, -0.5780, -0.1084, 1.0000], ])
CELL_COVARIANCE.loc[dict(group='IN_CC_L23')] = np.matrix([[1.0000, 0.1580, -0.5835, 0.4011, -0.0561, 0.0718, -0.2038, 0.2615, -0.2365,],
                              [0.1580, 1.0000, 0.0141, -0.1272, -0.4327, 0.1778, -0.0902, -0.0329, -0.3778,],
                              [-0.5835, 0.0141, 1.0000, -0.6295, -0.2949, -0.2008, 0.3164, -0.2615, -0.0536,],
                              [0.4011, -0.1272, -0.6295, 1.0000, 0.6960, -0.2587, -0.0988, 0.6113, 0.5636,],
                              [-0.0561, -0.4327, -0.2949, 0.6960, 1.0000, -0.3370, 0.2042, 0.3959, 0.8581,],
                              [0.0718, 0.1778, -0.2008, -0.2587, -0.3370, 1.0000, -0.0634, -0.5202, -0.3829,],
                              [-0.2038, -0.0902, 0.3164, -0.0988, 0.2042, -0.0634, 1.0000, 0.0559, 0.3322,],
                              [0.2615, -0.0329, -0.2615, 0.6113, 0.3959, -0.5202, 0.0559, 1.0000, 0.3210,],
                              [-0.2365, -0.3778, -0.0536, 0.5636, 0.8581, -0.3829, 0.3322, 0.3210, 1.0000,], ])
CELL_COVARIANCE.loc[dict(group='IN_F_L23')] = np.matrix([[1.0000, -0.1586, 0.1817, -0.0195, -0.0884, 0.0282, 0.0560, -0.1369, 0.0099,],
                             [-0.1586, 1.0000, 0.0440, 0.1013, -0.2510, -0.0046, -0.1105, 0.0738, -0.1152,],
                             [0.1817, 0.0440, 1.0000, -0.5118, 0.0414, 0.2570, 0.0932, 0.0961, 0.4938,],
                             [-0.0195, 0.1013, -0.5118, 1.0000, 0.0480, -0.1155, -0.2463, -0.0754, 0.0204,],
                             [-0.0884, -0.2510, 0.0414, 0.0480, 1.0000, 0.2577, -0.0581, 0.3152, 0.3151,],
                             [0.0282, -0.0046, 0.2570, -0.1155, 0.2577, 1.0000, -0.1598, 0.4397, 0.1107,],
                             [0.0560, -0.1105, 0.0932, -0.2463, -0.0581, -0.1598, 1.0000, -0.4617, 0.1872,],
                             [-0.1369, 0.0738, 0.0961, -0.0754, 0.3152, 0.4397, -0.4617, 1.0000, -0.0114,],
                             [0.0099, -0.1152, 0.4938, 0.0204, 0.3151, 0.1107, 0.1872, -0.0114, 1.0000,], ])
CELL_COVARIANCE.loc[dict(group='PC_L5')] = np.matrix([[1.0000, -0.2440, -0.2729, 0.2863, -0.0329, 0.2925, -0.0588, 0.3377, -0.1914,],
                           [-0.2440, 1.0000, 0.0874, -0.1523, -0.2565, -0.1605, 0.0874, -0.2895, -0.2125,],
                           [-0.2729, 0.0874, 1.0000, -0.6332, 0.2012, -0.0578, 0.0283, -0.1100, 0.3013,],
                           [0.2863, -0.1523, -0.6332, 1.0000, 0.3140, 0.2152, -0.1084, 0.4114, 0.1732,],
                           [-0.0329, -0.2565, 0.2012, 0.3140, 1.0000, 0.3184, -0.1923, 0.3761, 0.8433,],
                           [0.2925, -0.1605, -0.0578, 0.2152, 0.3184, 1.0000, 0.1246, 0.4736, 0.2078,],
                           [-0.0588, 0.0874, 0.0283, -0.1084, -0.1923, 0.1246, 1.0000, 0.0752, -0.1578,],
                           [0.3377, -0.2895, -0.1100, 0.4114, 0.3761, 0.4736, 0.0752, 1.0000, 0.2114,],
                           [-0.1914, -0.2125, 0.3013, 0.1732, 0.8433, 0.2078, -0.1578, 0.2114, 1.0000,], ])
CELL_COVARIANCE.loc[dict(group='IN_L_L5')] = np.matrix([[1.0000, -0.2894, 0.0381, 0.0664, -0.2418, 0.2253, 0.2822, -0.2919, 0.0581,],
                             [-0.2894, 1.0000, -0.2259, 0.4265, 0.1859, -0.6307, -0.0140, 0.4944, 0.2495,],
                             [0.0381, -0.2259, 1.0000, -0.2855, 0.0724, 0.1199, -0.1487, -0.3773, 0.1881,],
                             [0.0664, 0.4265, -0.2855, 1.0000, 0.2208, -0.3752, 0.0660, 0.3415, 0.7289,],
                             [-0.2418, 0.1859, 0.0724, 0.2208, 1.0000, 0.1412, -0.2931, 0.1993, 0.4609,],
                             [0.2253, -0.6307, 0.1199, -0.3752, 0.1412, 1.0000, -0.2855, -0.2046, -0.1974,],
                             [0.2822, -0.0140, -0.1487, 0.0660, -0.2931, -0.2855, 1.0000, -0.1172, -0.0851,],
                             [-0.2919, 0.4944, -0.3773, 0.3415, 0.1993, -0.2046, -0.1172, 1.0000, 0.0530,],
                             [0.0581, 0.2495, 0.1881, 0.7289, 0.4609, -0.1974, -0.0851, 0.0530, 1.0000,], ])
CELL_COVARIANCE.loc[dict(group='IN_L_d_L5')] = np.matrix([[1.0000, -0.2894, 0.0381, 0.0664, -0.2418, 0.2253, 0.2822, -0.2919, 0.0581,],
                               [-0.2894, 1.0000, -0.2259, 0.4265, 0.1859, -0.6307, -0.0140, 0.4944, 0.2495,],
                               [0.0381, -0.2259, 1.0000, -0.2855, 0.0724, 0.1199, -0.1487, -0.3773, 0.1881,],
                               [0.0664, 0.4265, -0.2855, 1.0000, 0.2208, -0.3752, 0.0660, 0.3415, 0.7289,],
                               [-0.2418, 0.1859, 0.0724, 0.2208, 1.0000, 0.1412, -0.2931, 0.1993, 0.4609,],
                               [0.2253, -0.6307, 0.1199, -0.3752, 0.1412, 1.0000, -0.2855, -0.2046, -0.1974,],
                               [0.2822, -0.0140, -0.1487, 0.0660, -0.2931, -0.2855, 1.0000, -0.1172, -0.0851,],
                               [-0.2919, 0.4944, -0.3773, 0.3415, 0.1993, -0.2046, -0.1172, 1.0000, 0.0530,],
                               [0.0581, 0.2495, 0.1881, 0.7289, 0.4609, -0.1974, -0.0851, 0.0530, 1.0000,], ])
CELL_COVARIANCE.loc[dict(group='IN_CL_L5')] = np.matrix([[1.0000, -0.2394, -0.6001, 0.3114, -0.2367, 0.5856, 0.2077, 0.0171, -0.4079,],
                              [-0.2394, 1.0000, -0.1764, 0.4675, 0.1810, -0.4942, -0.4389, 0.6950, 0.0811,],
                              [-0.6001, -0.1764, 1.0000, -0.6002, 0.2170, -0.0922, 0.2129, -0.3566, 0.4204,],
                              [0.3114, 0.4675, -0.6002, 1.0000, 0.2597, -0.1039, -0.5507, 0.7230, 0.0775,],
                              [-0.2367, 0.1810, 0.2170, 0.2597, 1.0000, 0.2159, -0.7123, 0.0193, 0.8494,],
                              [0.5856, -0.4942, -0.0922, -0.1039, 0.2159, 1.0000, 0.0587, -0.4724, 0.0957,],
                              [0.2077, -0.4389, 0.2129, -0.5507, -0.7123, 0.0587, 1.0000, -0.3395, -0.5780,],
                              [0.0171, 0.6950, -0.3566, 0.7230, 0.0193, -0.4724, -0.3395, 1.0000, -0.1084,],
                              [-0.4079, 0.0811, 0.4204, 0.0775, 0.8494, 0.0957, -0.5780, -0.1084, 1.0000,], ])
CELL_COVARIANCE.loc[dict(group='IN_CL_AC_L5')] = np.matrix([[1.0000, -0.2394, -0.6001, 0.3114, -0.2367, 0.5856, 0.2077, 0.0171, -0.4079,],
                                 [-0.2394, 1.0000, -0.1764, 0.4675, 0.1810, -0.4942, -0.4389, 0.6950, 0.0811,],
                                 [-0.6001, -0.1764, 1.0000, -0.6002, 0.2170, -0.0922, 0.2129, -0.3566, 0.4204,],
                                 [0.3114, 0.4675, -0.6002, 1.0000, 0.2597, -0.1039, -0.5507, 0.7230, 0.0775,],
                                 [-0.2367, 0.1810, 0.2170, 0.2597, 1.0000, 0.2159, -0.7123, 0.0193, 0.8494,],
                                 [0.5856, -0.4942, -0.0922, -0.1039, 0.2159, 1.0000, 0.0587, -0.4724, 0.0957,],
                                 [0.2077, -0.4389, 0.2129, -0.5507, -0.7123, 0.0587, 1.0000, -0.3395, -0.5780,],
                                 [0.0171, 0.6950, -0.3566, 0.7230, 0.0193, -0.4724, -0.3395, 1.0000, -0.1084,],
                                 [-0.4079, 0.0811, 0.4204, 0.0775, 0.8494, 0.0957, -0.5780, -0.1084, 1.0000,], ])
CELL_COVARIANCE.loc[dict(group='IN_CC_L5')] = np.matrix([[1.0000, -0.2440, -0.2729, 0.2863, -0.0329, 0.2925, -0.0588, 0.3377, -0.1914,],
                              [-0.2440, 1.0000, 0.0874, -0.1523, -0.2565, -0.1605, 0.0874, -0.2895, -0.2125,],
                              [-0.2729, 0.0874, 1.0000, -0.6332, 0.2012, -0.0578, 0.0283, -0.1100, 0.3013,],
                              [0.2863, -0.1523, -0.6332, 1.0000, 0.3140, 0.2152, -0.1084, 0.4114, 0.1732,],
                              [-0.0329, -0.2565, 0.2012, 0.3140, 1.0000, 0.3184, -0.1923, 0.3761, 0.8433,],
                              [0.2925, -0.1605, -0.0578, 0.2152, 0.3184, 1.0000, 0.1246, 0.4736, 0.2078,],
                              [-0.0588, 0.0874, 0.0283, -0.1084, -0.1923, 0.1246, 1.0000, 0.0752, -0.1578,],
                              [0.3377, -0.2895, -0.1100, 0.4114, 0.3761, 0.4736, 0.0752, 1.0000, 0.2114,],
                              [-0.1914, -0.2125, 0.3013, 0.1732, 0.8433, 0.2078, -0.1578, 0.2114, 1.0000,], ])
CELL_COVARIANCE.loc[dict(group='IN_F_L5')] = np.matrix([[1.0000, -0.1586, 0.1817, -0.0195, -0.0884, 0.0282, 0.0560, -0.1369, 0.0099,],
                             [-0.1586, 1.0000, 0.0440, 0.1013, -0.2510, -0.0046, -0.1105, 0.0738, -0.1152,],
                             [0.1817, 0.0440, 1.0000, -0.5118, 0.0414, 0.2570, 0.0932, 0.0961, 0.4938,],
                             [-0.0195, 0.1013, -0.5118, 1.0000, 0.0480, -0.1155, -0.2463, -0.0754, 0.0204,],
                             [-0.0884, -0.2510, 0.0414, 0.0480, 1.0000, 0.2577, -0.0581, 0.3152, 0.3151,],
                             [0.0282, -0.0046, 0.2570, -0.1155, 0.2577, 1.0000, -0.1598, 0.4397, 0.1107,],
                             [0.0560, -0.1105, 0.0932, -0.2463, -0.0581, -0.1598, 1.0000, -0.4617, 0.1872,],
                             [-0.1369, 0.0738, 0.0961, -0.0754, 0.3152, 0.4397, -0.4617, 1.0000, -0.0114,],
                             [0.0099, -0.1152, 0.4938, 0.0204, 0.3151, 0.1107, 0.1872, -0.0114, 1.0000,], ])                         
                  
        
K_TRANSF = np.zeros((len(MEMB_PARAM_NAMES), len(GROUP_NAMES)))
K_TRANSF = xr.DataArray(K_TRANSF, coords=[MEMB_PARAM_NAMES,GROUP_NAMES], dims=['param', 'group'], name='Memb_param K')

K_TRANSF[:,:]= np.matrix([[0.3700, 0.2200, 0.2200, 0.0000, 0.0000, 0.3700, 0.0000, 0.2300, 0.2200, 0.2200, 0.0000, 0.0000, 0.2300, 0.0000],
               [0.0000, 0.0200, 0.0200, 0.0000, 0.0000, 0.0000, 0.0100, 0.0100, 0.0200, 0.0200, 0.0000, 0.0000, 0.0100, 0.0100,],
               [0.0000, 0.3600, 0.3600, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3600, 0.3600, 0.0000, 0.0000, 0.0000, 0.0000,],
               [0.0100, 0.0000, 0.0000, 0.0000, 0.0000, 0.0100, 0.0000, 0.1300, 0.0000, 0.0000, 0.0000, 0.0000, 0.1300, 0.0000,],
               [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0200, 0.0000, 0.0000, 0.0000, 0.0000, 0.0200, 0.0000,],
               [0.0000, 0.0200, 0.0200, 0.0200, 0.0200, 0.0000, 0.0000, 0.0000, 0.0200, 0.0200, 0.0200, 0.0200, 0.0000, 0.0000,],
               [0.0100, 0.3600, 0.3600, 0.0000, 0.0000, 0.0100, 0.0000, 0.0000, 0.3600, 0.3600, 0.0000, 0.0000, 0.0000, 0.0000,],
               [0.0000, 0.0000, 0.0000, 0.1200, 0.1200, 0.0000, 0.2600, 0.0000, 0.0000, 0.0000, 0.1200, 0.1200, 0.0000, 0.2600,],
               [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,], ])



CELL_STD = np.zeros((len(MEMB_PARAM_NAMES), len(GROUP_NAMES)))
CELL_STD = xr.DataArray(CELL_STD, coords=[MEMB_PARAM_NAMES, GROUP_NAMES], dims=['param', 'group'], name='Memb_param std')

CELL_STD.loc[dict(param='C')] = [0.4296, 0.0754, 0.0754, 0.3283, 0.3283, 0.4296, 0.5029,
                                    0.1472, 0.0754, 0.0754, 0.3283, 0.3283, 0.1472, 0.5029 ] 
CELL_STD.loc[dict(param='g_L')] = [ 0.3558, 0.0046, 0.0046, 0.1844, 0.1844, 0.3558, 0.0022,
                                       0.0030, 0.0046, 0.0046, 0.1844, 0.1844, 0.0030, 0.0022 ]
CELL_STD.loc[dict(param='E_L')] = [ 0.3644, 0.3813, 0.3813, 0.3630, 0.3630, 0.3644, 0.3359,
                                       0.2846, 0.3813, 0.3813, 0.3630, 0.3630, 0.2846, 0.3359 ]
CELL_STD.loc[dict(param='delta_T')] = [ 0.0048, 0.7107, 0.7107, 0.3568, 0.3568, 0.0048, 0.7395,
                                           0.0554, 0.7107, 0.7107, 0.3568, 0.3568, 0.0554, 0.7395 ]
CELL_STD.loc[dict(param='V_up')] = [ 0.5259, 0.5033, 0.5033, 0.4372, 0.4372, 0.5259, 0.3035, 
                                        0.0062, 0.5033, 0.5033, 0.4372, 0.4372, 0.0062, 0.3035 ]
CELL_STD.loc[dict(param='tau_w')] = [ 0.4947, 0.0052, 0.0052, 0.0170, 0.0170, 0.4947, 0.3186, 
                                         0.6356, 0.0052, 0.0052, 0.0170, 0.0170, 0.6356, 0.3186 ]
CELL_STD.loc[dict(param='b')] = [ 0.0113, 1.9269, 1.9269, 1.4550, 1.4550, 0.0113, 1.0353, 
                                     1.3712, 1.9269, 1.9269, 1.4550, 1.4550, 1.3712, 1.0353 ]
CELL_STD.loc[dict(param='V_r')] = [ 0.6104, 0.4817, 0.4817, 0.1504, 0.1504, 0.6104, 0.1813, 
                                       0.3497, 0.4817, 0.4817, 0.1504, 0.1504, 0.3497, 0.1813 ]
CELL_STD.loc[dict(param='V_T')] = [ 0.4608, 0.4385, 0.4385, 0.4311, 0.4311, 0.4608, 0.3632, 
                                       0.2857, 0.4385, 0.4385, 0.4311, 0.4311, 0.2857, 0.3632 ]

CELL_MINIMUM = np.zeros((len(MEMB_PARAM_NAMES), len(GROUP_NAMES)))
CELL_MINIMUM = xr.DataArray(CELL_MINIMUM, coords=[MEMB_PARAM_NAMES, GROUP_NAMES], dims=['param', 'group'], name='Memb_param min')

CELL_MINIMUM.loc[dict(param='C')] = [61.4187,42.1156,42.1156,51.8447,51.8447,61.4187,32.3194,
                                        110.7272,42.1156,42.1156,51.8447,51.8447,110.7272,32.3194,]
CELL_MINIMUM.loc[dict(param='g_L')] = [3.2940,3.6802,3.6802,2.9852,2.9852,3.2940,2.1462,
                                          3.4510,3.6802,3.6802,2.9852,2.9852,3.4510,2.1462,]
CELL_MINIMUM.loc[dict(param='E_L')] = [-104.9627,-96.9345,-96.9345,-98.8335,-98.8335,-104.9627,-102.3895,
                                          -101.5624,-96.9345,-96.9345,-98.8335,-98.8335,-101.5624,-102.3895,]
CELL_MINIMUM.loc[dict(param='delta_T')] = [10.5568,2.1840,2.1840,11.0503,11.0503,10.5568,1.8285,
                                              12.7969,2.1840,2.1840,11.0503,11.0503,12.7969,1.8285,]
CELL_MINIMUM.loc[dict(param='V_up')] = [-62.5083,-60.6745,-60.6745,-65.4193,-65.4193,-62.5083,-42.8895,
                                           -66.1510,-60.6745,-60.6745,-65.4193,-65.4193,-66.1510,-42.8895,]
CELL_MINIMUM.loc[dict(param='tau_w')] = [54.0018,10.2826,10.2826,12.2898,12.2898,54.0018,20.0311,
                                            33.1367,10.2826,10.2826,12.2898,12.2898,33.1367,20.0311,]
CELL_MINIMUM.loc[dict(param='b')] = [1.2406,1.0000,1.0000,1.0000,1.0000,1.2406,1.0000,
                                        1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,]
CELL_MINIMUM.loc[dict(param='V_r')] = [-219.2039,-128.4559,-128.4559,-271.9846,-271.9846,-219.2039,-105.1880,
                                          -124.5158,-128.4559,-128.4559,-271.9846,-271.9846,-124.5158,-105.1880,]
CELL_MINIMUM.loc[dict(param='V_T')] = [-63.2375,-85.2096,-85.2096,-70.3537,-70.3537,-63.2375,-53.3897,
                                          -69.5922,-85.2096,-85.2096,-70.3537,-70.3537,-69.5922,-53.3897,]


CELL_MAXIMUM = np.zeros((len(MEMB_PARAM_NAMES), len(GROUP_NAMES)))
CELL_MAXIMUM = xr.DataArray(CELL_MAXIMUM, coords=[MEMB_PARAM_NAMES, GROUP_NAMES], dims=['param', 'group'], name='Memb_param max')

CELL_MAXIMUM.loc[dict(param='C')] = [337.9765,94.6939,94.6939,126.2367,126.2367,337.9765,201.3221,
                                        617.2776,94.6939,94.6939,126.2367,126.2367,617.2776,201.3221,]
CELL_MAXIMUM.loc[dict(param='g_L')] = [10.8106,8.6130,8.6130,5.6192,5.6192,10.8106,5.3460, 
                                          15.6329,8.6130,8.6130,5.6192,5.6192,15.6329,5.3460,]
CELL_MAXIMUM.loc[dict(param='E_L')] = [-76.8526,-71.7548,-71.7548,-75.7868,-75.7868,-76.8526,-59.6898, 
                                          -66.4770,-71.7548,-71.7548,-75.7868,-75.7868,-66.4770,-59.6898,]
CELL_MAXIMUM.loc[dict(param='delta_T')] = [45.3814,40.4333,40.4333,31.3533,31.3533,45.3814,47.6214, 
                                              43.5882,40.4333,40.4333,31.3533,31.3533,43.5882,47.6214,]
CELL_MAXIMUM.loc[dict(param='V_up')] = [-30.0577,-36.5929,-36.5929,-45.6445,-45.6445,-30.0577,-30.7977, 
                                           -25.2891,-36.5929,-36.5929,-45.6445,-45.6445,-25.2891,-30.7977,]
CELL_MAXIMUM.loc[dict(param='tau_w')] = [232.8699,21.9964,21.9964,120.5043,120.5043,232.8699,102.4180, 
                                            909.5520,21.9964,21.9964,120.5043,120.5043,909.5520,102.4180,]
CELL_MAXIMUM.loc[dict(param='b')] = [40.2930,196.7634,196.7634,71.0958,71.0958,40.2930,54.2781, 
                                        325.7906,196.7634,196.7634,71.0958,71.0958,325.7906,54.2781,]
CELL_MAXIMUM.loc[dict(param='V_r')] = [-45.0393,-56.5047,-56.5047,-56.8682,-56.8682,-45.0393,-35.7409, 
                                          -35.1145,-56.5047,-56.5047,-56.8682,-56.8682,-35.1145,-35.7409,]
CELL_MAXIMUM.loc[dict(param='V_T')] = [-36.8701,-39.1085,-39.1085,-49.0974,-49.0974,-36.8701,-20.6720, 
                                          -27.8669,-39.1085,-39.1085,-49.0974,-49.0974,-27.8669,-20.6720,]

TAU_MIN = xr.DataArray(np.asarray([10.3876,7.3511,7.3511,9.2264,9.2264,10.3876,5.8527,
		16.7015,7.3511,7.3511,9.2264,9.2264,16.7015,5.8527,]),coords=[GROUP_NAMES,], dims='group', name='Memb_param tau_m min')   

TAU_MAX =  xr.DataArray(np.asarray([42.7304,15.9128,15.9128,25.9839,25.9839,42.7304,48.7992,
		67.7062,15.9128,15.9128,25.9839,25.9839,67.7062,48.7992,]),coords=[GROUP_NAMES,], dims='group', name='Memb_param tau_m max')


INTERSTRIPE_SETS ={}

INTERSTRIPE_SETS['A']={
    'pair': ['PC_L23','PC_L23'],
    'coefficient_0': 0.909,
    'coefficient_1': 1.4587,
    'connections': [-4, -2, 2, 4],
    }

INTERSTRIPE_SETS['B']={
    'pair': ['PC_L5','IN_CC_L5'],
    'coefficient_0': 0.909,
    'coefficient_1': 1.0375,
    'connections':[-1, 1],
    }

INTERSTRIPE_SETS['C']={
    'pair':  ['PC_L5', 'IN_F_L5'],
    'coefficient_0': 0.909,
    'coefficient_1': 1.0375,
    'connections': [-1, 1],
    }

PCON = np.zeros((len(GROUP_NAMES), len(GROUP_NAMES)))
PCON = xr.DataArray(PCON, coords=[GROUP_NAMES, GROUP_NAMES], dims=['target', 'source'], name='pCon')


PCON.loc[dict(target='PC_L23', source='PC_L23')] = 0.1393
PCON.loc[dict(target='PC_L23', source='PC_L5')] = 0.0449
PCON.loc[dict(target='PC_L5', source='PC_L23')] = 0.2333
PCON.loc[dict(target='PC_L5', source='PC_L5')] = 0.0806

PCON.loc[dict(target=['IN_L_L23', 'IN_L_d_L23'], source='PC_L23')] = 0.3247
PCON.loc[dict(target=['IN_L_L23', 'IN_L_d_L23'], source='PC_L5')] = 0.1875
PCON.loc[dict(target=['IN_L_L5', 'IN_L_d_L5'], source='PC_L23')] = 0.0870
PCON.loc[dict(target=['IN_L_L5', 'IN_L_d_L5'], source='PC_L5')] = 0.3331  
PCON.loc[dict(target=['IN_CL_L23', 'IN_CL_AC_L23'], source='PC_L23')] = 0.1594
PCON.loc[dict(target=['IN_CL_L23', 'IN_CL_AC_L23'], source='PC_L5')] = 0.0920
PCON.loc[dict(target=['IN_CL_L5', 'IN_CL_AC_L5'], source='PC_L23')] = 0.0800
PCON.loc[dict(target=['IN_CL_L5', 'IN_CL_AC_L5'], source='PC_L5')] = 0.0800   

PCON.loc[dict(target='IN_CC_L23', source='PC_L23')] = 0.3247
PCON.loc[dict(target='IN_CC_L23', source='PC_L5')] = 0.1875
PCON.loc[dict(target='IN_CC_L5', source='PC_L23')] = 0.0870
PCON.loc[dict(target='IN_CC_L5', source='PC_L5')] = 0.3331
PCON.loc[dict(target='IN_F_L23', source='PC_L23')] = 0.2900
PCON.loc[dict(target='IN_F_L23', source='PC_L5')] = 0.1674
PCON.loc[dict(target='IN_F_L5', source='PC_L23')] = 0.1500
PCON.loc[dict(target='IN_F_L5', source='PC_L5')] = 0.3619

PCON.loc[dict(target='PC_L23', source=['IN_L_L23', 'IN_L_d_L23'])] = 0.4586
PCON.loc[dict(target='PC_L23', source=['IN_L_L5', 'IN_L_d_L5'])] = 0.0991
PCON.loc[dict(target='PC_L5', source=['IN_L_L23', 'IN_L_d_L23'])] = 0.2130
PCON.loc[dict(target='PC_L5', source=['IN_L_L5', 'IN_L_d_L5'])] = 0.7006
PCON.loc[dict(target='PC_L23', source=['IN_CL_L23', 'IN_CL_AC_L23'])] = 0.4164
PCON.loc[dict(target='PC_L23', source=['IN_CL_L5', 'IN_CL_AC_L5'])] = 0.0321
PCON.loc[dict(target='PC_L5', source=['IN_CL_L23', 'IN_CL_AC_L23'])] = 0.1934
PCON.loc[dict(target='PC_L5', source=['IN_CL_L5', 'IN_CL_AC_L5'])] = 0.2271
PCON.loc[dict(target='PC_L23', source='IN_CC_L23')] = 0.4586
PCON.loc[dict(target='PC_L23', source='IN_CC_L5')] = 0.0991
PCON.loc[dict(target='PC_L5', source='IN_CC_L23')] = 0.2130
PCON.loc[dict(target='PC_L5', source='IN_CC_L23')] = 0.7006
PCON.loc[dict(target='PC_L23', source='IN_F_L23')] = 0.6765
PCON.loc[dict(target='PC_L23', source='IN_F_L5')] = 0.1287
PCON.loc[dict(target='PC_L5', source='IN_F_L23')] = 0.3142
PCON.loc[dict(target='PC_L5', source='IN_F_L5')] = 0.9096

 

PCON.loc[dict(target=GROUPS_SETS['IN_L23'], source=GROUPS_SETS['IN_L23'])]   = 0.25
PCON.loc[dict(target=GROUPS_SETS['IN_L5'], source=GROUPS_SETS['IN_L5'])]    = 0.60


CLUSTER_FLAG = np.zeros((14,14))
CLUSTER_FLAG = xr.DataArray(CLUSTER_FLAG, coords=[GROUP_NAMES,GROUP_NAMES], dims=['target', 'source'], name='Clustering flag')                      

CLUSTER_FLAG.loc[dict(target='PC_L23', source='PC_L23')] = 1
CLUSTER_FLAG.loc[dict(target='PC_L5', source='PC_L5')] = 1


STSP_PARAMS = {'U': dict(unit=1, value='U'), 'tau_rec': dict(unit='ms', value='tau_rec'), 'tau_fac': dict(unit='ms', value='tau_rec')}
STSP_VARS = {'u_temp': dict(unit=1, value='U'), 'u': dict(unit=1, value='U'),
             'R_temp': dict(unit=1, value=1),'R': dict(unit=1, value=1) ,
             'a_syn':dict(unit=1, value='U')}
STSP_DECL_VARS = STSP_PARAMS | STSP_VARS

STSP_PARAM_NAMES = list(STSP_PARAMS.keys())
STSP_VAR_NAMES = list(STSP_VARS.keys())
STSP_DECL_VAR_NAMES = list(STSP_DECL_VARS.keys())

STSP_PARAM_DISTR= ['{}_{}'.format(var, stat) for stat in ['mean', 'std'] for var in STSP_PARAM_NAMES]


SYN_PARAMS = {'gmax': dict(unit='nS', value='gmax'),
              'pfail': dict(unit=1, value='pfail')}
SYN_PARAM_NAMES = list(SYN_PARAMS.keys())

SYN_AUX = {'failure': dict(unit=1), 'gsyn_amp': dict(unit=1)}
SYN_AUX_NAMES = list(SYN_AUX.keys())

STSP_KINDS_DISTR_DICT = {
    'E1': np.array([0.28, 194, 507, 0.02, 18, 37]),
    'E2': np.array([0.25, 671, 17, 0.02, 17, 5]),
    'E3': np.array([0.29, 329, 326, 0.03, 53, 66]),
    'I1': np.array([0.16, 45, 376,0.10, 21, 253]),
    'I2': np.array([0.25, 706, 21, 0.13, 405, 9]),
    'I3': np.array([0.32, 144, 62, 0.14, 80, 31]),
    }


STPS_KINDS = xr.DataArray(np.array(list(STSP_KINDS_DISTR_DICT.values())), coords=[list(STSP_KINDS_DISTR_DICT.keys()), STSP_PARAM_DISTR], dims=['kind', 'param'], name='STSP param sets')


STSP_SETS_DICT={
    'A': np.asarray([0.45, 0.38, 0.17, np.NaN, np.NaN, np.NaN]),
    'B': np.asarray([1, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]),
    'C': np.asarray([np.NaN, 1, np.NaN, np.NaN, np.NaN, np.NaN]),
    'D': np.asarray([np.NaN, np.NaN, np.NaN, 0.25, 0.5, 0.25]),
    'E': np.asarray([np.NaN, np.NaN, np.NaN, np.NaN, 1, np.NaN]),
    'F': np.asarray([np.NaN, np.NaN, np.NaN, 0.29, 0.58, 0.13])
    }

STSP_SETS = xr.DataArray(np.array(list(STSP_SETS_DICT.values())), coords=[list(STSP_SETS_DICT.keys()), list(STSP_KINDS_DISTR_DICT.keys())], dims=['set', 'kind'], name='STSP set distribution')

STSP_SET_DISTR = np.asarray([['' for i in range(len(GROUP_NAMES))]\
                         for j in range(len(GROUP_NAMES))], dtype='U16')

STSP_SET_DISTR = xr.DataArray(STSP_SET_DISTR, coords=[GROUP_NAMES, GROUP_NAMES], dims=['target', 'source'], name='STSP distribution groups')
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['PC'], source=GROUPS_SETS['PC'])] = 'A'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['IN_L'], source=GROUPS_SETS['PC'])] = 'B'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['IN_L_d'], source=GROUPS_SETS['PC'])] = 'C'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['IN_CL'], source=GROUPS_SETS['PC'])] = 'C'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['IN_CL_AC'], source=GROUPS_SETS['PC'])] = 'B'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['IN_CC'], source=GROUPS_SETS['PC'])] = 'C'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['IN_F'], source=GROUPS_SETS['PC'])] = 'B'

STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['PC'], source=GROUPS_SETS['IN_L'])] = 'D'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['PC'], source=GROUPS_SETS['IN_L_d'])] = 'E'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['PC'], source=GROUPS_SETS['IN_CL'])] = 'E'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['PC'], source=GROUPS_SETS['IN_CL_AC'])] = 'E'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['PC'], source=GROUPS_SETS['IN_CC'])] = 'E'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['PC'], source=GROUPS_SETS['IN_F'])] = 'E'

STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['IN_L23'], source=GROUPS_SETS['IN_L23'])] = 'F'
STSP_SET_DISTR.loc[dict(target=GROUPS_SETS['IN_L5'], source=GROUPS_SETS['IN_L5'])] = 'F'



SYN_DELAY = np.zeros((len(GROUP_NAMES), len(GROUP_NAMES)))
SYN_DELAY = xr.DataArray(SYN_DELAY, coords=[GROUP_NAMES, GROUP_NAMES], dims=['target', 'source'], name='Syn delay')

SYN_DELAY.loc[dict(target='PC_L23', source='PC_L23')] = 1.5465
SYN_DELAY.loc[dict(target='PC_L23', source='PC_L5')] = 2.7533
SYN_DELAY.loc[dict(target='PC_L5', source='PC_L23')] = 1.9085
SYN_DELAY.loc[dict(target='PC_L5', source='PC_L5')] = 1.5667

SYN_DELAY.loc[dict(target='PC_L23', source=GROUPS_SETS['IN_L23'])] = 1.2491
SYN_DELAY.loc[dict(target='PC_L23', source=GROUPS_SETS['IN_L5'])] = 1.4411
SYN_DELAY.loc[dict(target='PC_L5', source=GROUPS_SETS['IN_L23'])] = 1.5415  
SYN_DELAY.loc[dict(target='PC_L5', source=GROUPS_SETS['IN_L5'])]  = 0.82

SYN_DELAY.loc[dict(target=GROUPS_SETS['IN_L23'], source='PC_L23')] = 0.9581
SYN_DELAY.loc[dict(target=GROUPS_SETS['IN_L23'], source='PC_L5')]  = 1.0544
SYN_DELAY.loc[dict(target=GROUPS_SETS['IN_L5'], source='PC_L23')]  = 1.1825
SYN_DELAY.loc[dict(target=GROUPS_SETS['IN_L5'], source='PC_L5')] = 0.6
SYN_DELAY.loc[dict(target=GROUPS_SETS['IN_L23'], source=GROUPS_SETS['IN_L23'])]   = 1.1
SYN_DELAY.loc[dict(target=GROUPS_SETS['IN_L5'], source=GROUPS_SETS['IN_L5'])]    = 1.1

SYN_GMAX = np.zeros((len(GROUP_NAMES), len(GROUP_NAMES)))
SYN_GMAX = xr.DataArray(SYN_GMAX, coords=[GROUP_NAMES, GROUP_NAMES], dims=['target', 'source'], name='Syn gmax')
  
SYN_GMAX.loc[dict(target='PC_L23', source=GROUPS_SETS['ALL'])] = [
    0.8405, 2.2615, 2.2615, 0.18, 0.18, 2.2615, 1.8218,
    0.8378, 0.2497, 0.2497, 0.0556, 0.0556, 0.2497, 0.2285,
    ]

SYN_GMAX.loc[dict(target='PC_L5', source=GROUPS_SETS['ALL'])]  = [
    0.9533, 1.0503, 1.0503, 0.0836, 0.0836, 1.0503, 0.8461, 
    0.8818, 1.7644, 1.7644, 0.3932, 0.3932, 1.7644, 1.6146,
    ]

SYN_GMAX.loc[dict(target=GROUPS_SETS['IN_L23'], source='PC_L23')] = [
    1.3403, 1.3403, 0.4710, 0.4710, 1.3403, 0.2500
    ]
SYN_GMAX.loc[dict(target=GROUPS_SETS['IN_L5'], source='PC_L23')]  = [
    1.5201, 1.5201, 0.5342, 0.5342, 1.5201, 0.2835
    ]

SYN_GMAX.loc[dict(target=GROUPS_SETS['IN_L23'], source='PC_L5')]  = [
    0.7738, 0.7738, 0.2719, 0.2719, 0.7738, 0.1443
    ]


SYN_GMAX.loc[dict(target=GROUPS_SETS['IN_L5'], source='PC_L5')] = [
    1.7431, 1.7431, 0.88, 0.88, 1.7431, 0.28
    ]

SYN_GMAX.loc[dict(target=GROUPS_SETS['IN_L23'], source=GROUPS_SETS['IN_L23'])]   = 1.35    
SYN_GMAX.loc[dict(target=GROUPS_SETS['IN_L5'], source=GROUPS_SETS['IN_L5'])]    = 1.35

SYN_SOURCE_GMAX_CORRECTION = {
    'PC':np.array([[
        1.0569, 0.5875, 0.6587, 0.7567, 0.6728, 0.9899, 0.6294,
        1.6596, 0.5941, 0.6661, 0.7647, 0.6799, 1.5818, 0.6360
        ]]).transpose(),
    'IN':np.array([[
        2.3859, 1.6277, 1.6277, 1.6671, 1.6671, 2.3142, 1.4363, 
        3.5816, 1.6277, 1.6277, 1.6671, 1.6671, 3.4016, 1.4363
        ]]).transpose()
}

for kind in ['PC', 'IN']: 
    SYN_GMAX.loc[dict(target=GROUPS_SETS['ALL'], source=GROUPS_SETS[kind])] = SYN_GMAX.loc[dict(target=GROUPS_SETS['ALL'], source=GROUPS_SETS[kind])]*repmat(SYN_SOURCE_GMAX_CORRECTION[kind], 1, len(GROUPS_SETS[kind]))

PFAIL=0.3
SYN_PFAIL = np.ones((len(GROUP_NAMES), len(GROUP_NAMES)))*PFAIL
SYN_PFAIL = xr.DataArray(SYN_PFAIL, coords=[GROUP_NAMES, GROUP_NAMES], dims=['target', 'source'], name='Syn pfail')

SYN_KINDS = np.array([[tp for tp in GROUP_KINDS.values()] for i in range(len(GROUP_NAMES))])
SYN_KINDS = xr.DataArray(SYN_KINDS, coords=[GROUP_NAMES, GROUP_NAMES], dims=['target', 'source'], name='Syn kinds')


AMPA = np.array([1.4, 10, 0, 0, 0, 0], dtype='float64') # in ms and mV
GABA = np.asarray([3, 40, -70, 0, 0, 0], dtype='float64') # in ms and mV
NMDA = np.array([4.3, 75, 0, 1, 0.0625, 0], dtype='float64') # in ms and mV

PARAM_CHANNELS = xr.DataArray(np.asarray([AMPA, GABA, NMDA]), coords=[list(CHANNEL_NAMES), CHANNEL_PARAMS], dims=['channel', 'param'])

GSYN_FACTOR_CHANNELS = xr.DataArray(np.asarray([[1],[1],[1.09]]), coords=[list(CHANNEL_NAMES), ['factor']], dims=['channel', 'param'])

GMAX_SIGMA = np.zeros((len(GROUP_NAMES), len(GROUP_NAMES)))   
GMAX_SIGMA = xr.DataArray(GMAX_SIGMA, coords=[GROUP_NAMES, GROUP_NAMES],dims=['target', 'source'])      

GMAX_SIGMA.loc[dict(target='PC_L23', source='PC_L23')] = 0.4695
GMAX_SIGMA.loc[dict(target='PC_L23', source='PC_L5')] = 0.1375
GMAX_SIGMA.loc[dict(target='PC_L5', source='PC_L23')] = 0.3530
GMAX_SIGMA.loc[dict(target='PC_L5', source='PC_L5')] = 0.9653

GMAX_SIGMA.loc[dict(target=['IN_L_L23', 'IN_L_d_L23'], source='PC_L23')] = 1.0855
GMAX_SIGMA.loc[dict(target=['IN_L_L23', 'IN_L_d_L23'], source='PC_L5')] = 0.6267
GMAX_SIGMA.loc[dict(target=['IN_L_L5', 'IN_L_d_L5'], source='PC_L23')] = 0.8588
GMAX_SIGMA.loc[dict(target=['IN_L_L5', 'IN_L_d_L5'], source='PC_L5')] = 1.1194

GMAX_SIGMA.loc[dict(target=['IN_CL_L23', 'IN_CL_AC_L23'], source='PC_L23')] = 0.1999
GMAX_SIGMA.loc[dict(target=['IN_CL_L23', 'IN_CL_AC_L23'], source='PC_L5')] = 0.1154
GMAX_SIGMA.loc[dict(target=['IN_CL_L5', 'IN_CL_AC_L5'], source='PC_L23')] = 0.1581
GMAX_SIGMA.loc[dict(target=['IN_CL_L5', 'IN_CL_AC_L5'], source='PC_L5')] = 0.7033

GMAX_SIGMA.loc[dict(target='IN_CC_L23', source='PC_L23')] = 1.0855
GMAX_SIGMA.loc[dict(target='IN_CC_L23', source='PC_L5')] = 0.6267
GMAX_SIGMA.loc[dict(target='IN_CC_L5', source='PC_L23')] = 0.8588
GMAX_SIGMA.loc[dict(target='IN_CC_L5', source='PC_L5')] = 1.1194

GMAX_SIGMA.loc[dict(target='IN_F_L23', source='PC_L23')] = 0.2000
GMAX_SIGMA.loc[dict(target='IN_F_L23', source='PC_L5')] = 0.1155
GMAX_SIGMA.loc[dict(target='IN_F_L5', source='PC_L23')] = 0.1582
GMAX_SIGMA.loc[dict(target='IN_F_L5', source='PC_L5')] = 0.3000

GMAX_SIGMA.loc[dict(target='PC_L23', source=['IN_L_L23', 'IN_L_d_L23'])] = 1.9462
GMAX_SIGMA.loc[dict(target='PC_L23', source=['IN_L_L5', 'IN_L_d_L5'])] = 0.0362
GMAX_SIGMA.loc[dict(target='PC_L5', source=['IN_L_L23', 'IN_L_d_L23'])] = 0.9038
GMAX_SIGMA.loc[dict(target='PC_L5', source=['IN_L_L5', 'IN_L_d_L5'])] = 0.2557

GMAX_SIGMA.loc[dict(target='PC_L23', source=['IN_CL_L23', 'IN_CL_AC_L23'])] = 0.6634
GMAX_SIGMA.loc[dict(target='PC_L23', source=['IN_CL_L5', 'IN_CL_AC_L5'])] = 0.0093
GMAX_SIGMA.loc[dict(target='PC_L5', source=['IN_CL_L23', 'IN_CL_AC_L23'])] = 0.3081
GMAX_SIGMA.loc[dict(target='PC_L5', source=['IN_CL_L5', 'IN_CL_AC_L5'])] = 0.0655

GMAX_SIGMA.loc[dict(target='PC_L23', source='IN_CC_L23')] = 1.9462
GMAX_SIGMA.loc[dict(target='PC_L23', source='IN_CC_L5')] = 0.0362
GMAX_SIGMA.loc[dict(target='PC_L5', source='IN_CC_L23')] = 0.9038
GMAX_SIGMA.loc[dict(target='PC_L5', source='IN_CC_L23')] = 0.2557

GMAX_SIGMA.loc[dict(target='PC_L23', source='IN_F_L23')] = 3.6531
GMAX_SIGMA.loc[dict(target='PC_L23', source='IN_F_L5')] = 0.1828
GMAX_SIGMA.loc[dict(target='PC_L5', source='IN_F_L23')] = 1.6966
GMAX_SIGMA.loc[dict(target='PC_L5', source='IN_F_L5')] = 1.2919

GMAX_SIGMA.loc[dict(target=GROUPS_SETS['IN_L23'], source=GROUPS_SETS['IN_L23'])]   = 0.35
GMAX_SIGMA.loc[dict(target=GROUPS_SETS['IN_L5'], source=GROUPS_SETS['IN_L5'])]    = 0.35


GMAX_MIN = np.zeros((len(GROUP_NAMES), len(GROUP_NAMES)))
GMAX_MIN = xr.DataArray(GMAX_MIN, coords=[GROUP_NAMES, GROUP_NAMES],dims=['target', 'source'])

GMAX_MAX = np.zeros((len(GROUP_NAMES), len(GROUP_NAMES)))
GMAX_MAX = xr.DataArray(GMAX_MAX, coords=[GROUP_NAMES, GROUP_NAMES],dims=['target', 'source'])
GMAX_MAX[:,:] = 100    
         

DELAY_SIGMA = np.zeros((len(GROUP_NAMES), len(GROUP_NAMES)))   
DELAY_SIGMA = xr.DataArray(DELAY_SIGMA, coords=[GROUP_NAMES, GROUP_NAMES],dims=['target', 'source'])

DELAY_SIGMA.loc[dict(target='PC_L23', source='PC_L23')] = 0.3095
DELAY_SIGMA.loc[dict(target='PC_L23', source='PC_L5')] = 0.1825
DELAY_SIGMA.loc[dict(target='PC_L5', source='PC_L23')] = 0.1651
DELAY_SIGMA.loc[dict(target='PC_L5', source='PC_L5')] = 0.4350  

DELAY_SIGMA.loc[dict(target=GROUPS_SETS['IN_L23'], source='PC_L23')] = 0.2489
DELAY_SIGMA.loc[dict(target=GROUPS_SETS['IN_L23'], source='PC_L5')]  = 0.0839
DELAY_SIGMA.loc[dict(target=GROUPS_SETS['IN_L5'], source='PC_L23')]  = 0.1327
DELAY_SIGMA.loc[dict(target=GROUPS_SETS['IN_L5'], source='PC_L5')] = 0.2000


DELAY_SIGMA.loc[dict(target='PC_L23', source=GROUPS_SETS['IN_L23'])] = 0.1786
DELAY_SIGMA.loc[dict(target='PC_L23', source=GROUPS_SETS['IN_L5'])] = 0.0394
DELAY_SIGMA.loc[dict(target='PC_L5', source=GROUPS_SETS['IN_L23'])] = 0.0940
DELAY_SIGMA.loc[dict(target='PC_L5', source=GROUPS_SETS['IN_L5'])] = 0.0940

DELAY_SIGMA.loc[dict(target=GROUPS_SETS['IN_L23'], source=GROUPS_SETS['IN_L23'])]   = 0.4
DELAY_SIGMA.loc[dict(target=GROUPS_SETS['IN_L5'], source=GROUPS_SETS['IN_L5'])]    = 0.4

DELAY_MIN = np.zeros((len(GROUP_NAMES), len(GROUP_NAMES)))   
DELAY_MIN = xr.DataArray(DELAY_MIN, coords=[GROUP_NAMES, GROUP_NAMES],dims=['target', 'source'])

DELAY_MAX = np.zeros((len(GROUP_NAMES), len(GROUP_NAMES)))   
DELAY_MAX = xr.DataArray(DELAY_MAX, coords=[GROUP_NAMES, GROUP_NAMES],dims=['target', 'source'])

DELAY_MAX[:,:] = 2    


spiking_params = {}
for par in SYN_PARAM_NAMES:
    spiking_params[par] = dict()
spiking_params['pfail'] = SYN_PFAIL
spiking_params['gmax'] = dict(mean=SYN_GMAX, sigma=GMAX_SIGMA, min=GMAX_MIN, max=GMAX_MAX)

EQ_MEMB_PARAM=['{}: {}'.format(name, UNIT_MAIN_DICT[MEMB_PARAM_UNITS[name]['unit']]) for name in MEMB_PARAM_NAMES]
EQ_MEMB_PARAM = '\n'.join(EQ_MEMB_PARAM)
CHANNEL_PARAMS_KINDS = ['{}_{}'.format(param,kind) for param in CHANNEL_PARAMS for kind in CHANNEL_KINDS]

# _=['{} = {}*{}: {}'.format(MEMB_PARAM_NAMES[i], list(MEMB_PARAM_UNITS.values())[i][0]) for i in range(len(MEMB_PARAM_NAMES))]
EQ_CHANNEL_PARAM=['{} = {}*{}: {}'.format('{}_{}'.format(param,kind), float(PARAM_CHANNELS.loc[kind,param].values), *CHANNEL_PARAMS_UNITS[param]) for param in CHANNEL_PARAMS for kind in CHANNEL_NAMES]
EQ_CHANNEL_PARAM = '\n'.join(EQ_CHANNEL_PARAM)

EQ_AUX_VAR =\
"""I_ref: amp
last_spike: second"""

EQ_MG_FACTOR =\
'Mg_{0} = (1/(1 +  Mg_fac_{0} * exp(Mg_slope_{0} * (Mg_half_{0}*mV - V)/mV))):1'
EQ_MG_FACTOR = [EQ_MG_FACTOR.format(name) for name in CHANNEL_NAMES]
EQ_MG_FACTOR = '\n'.join(EQ_MG_FACTOR)

EQ_CHANNEL_CURRENT =\
'I_{0} = g_{0} * (E_{0} - V) * Mg_{0}: amp'
EQ_CHANNEL_CURRENT = [EQ_CHANNEL_CURRENT.format(name) for name in CHANNEL_NAMES]
EQ_CHANNEL_CURRENT = '\n'.join(EQ_CHANNEL_CURRENT)

EQ_CHANNEL_CONDUCTANCE =\
'g_{0} = g_{0}_off - g_{0}_on: siemens'
EQ_CHANNEL_CONDUCTANCE = [EQ_CHANNEL_CONDUCTANCE.format(name) for name in CHANNEL_NAMES]
EQ_CHANNEL_CONDUCTANCE = '\n'.join(EQ_CHANNEL_CONDUCTANCE)

EQ_CHANNEL_DGDT  =\
"""dg_{0}_off/dt = - (1/tau_off_{0}) * g_{0}_off: siemens
dg_{0}_on/dt = - (1/tau_on_{0}) * g_{0}_on: siemens"""
EQ_CHANNEL_DGDT = [EQ_CHANNEL_DGDT.format(name) for name in CHANNEL_NAMES]
EQ_CHANNEL_DGDT = '\n'.join(EQ_CHANNEL_DGDT)

EQ_CURRENT =\
"""I_DC: amp
I_AC = {}*pA: amp""" + '\n'
EQ_CURRENT +=\
'I_syn = '+' + '.join(['I_{0}'.format(name) for name in CHANNEL_NAMES])+': amp\n'
EQ_CURRENT+=\
"""I_inj = I_DC + I_AC: amp
I_tot =  I_syn + I_inj: amp"""


EQ_MEMB = \
"""I_exp = g_L * delta_T * exp((V - V_T)/delta_T): amp
w_V = I_tot + I_exp -g_L * (V - E_L): amp
dV = int(I_tot >= I_ref) * int(t - last_spike < 5 * ms) * (-g_L/C)*(V - V_r) + (1 - int(I_tot >= I_ref) * int(t - last_spike < 5 * ms)) * (I_tot + I_exp - g_L * (V - E_L) - w)/C: volt/second
dV/dt = dV: volt

D0 = (C/g_L) * w_V:  coulomb
dD0 = C *(exp((V - V_T)/delta_T)-1): farad
dw/dt = int(w > w_V - D0/tau_w) * int(w < w_V + D0/tau_w) * int(V <= V_T) * int(I_tot < I_ref) * -(g_L * (1 - exp((V - V_T)/delta_T)) + dD0/tau_w)*dV: amp"""


EQ_MEMBR_MODEL = '\n\n'.join([EQ_MEMB_PARAM, EQ_CHANNEL_PARAM,
                         EQ_MG_FACTOR,EQ_CHANNEL_CURRENT,EQ_CHANNEL_CONDUCTANCE,
                         EQ_CHANNEL_DGDT,  EQ_AUX_VAR, EQ_CURRENT, EQ_MEMB])

EQ_MEMBR_THRESHOLD = "V > V_up"
EQ_MEMBR_RESET = "V = V_r; w += b"
EQ_MEMBR_EVENT = {'w_crossing':{}}
EQ_MEMBR_EVENT['w_crossing']['condition'] = 'w > w_V - D0/tau_w and w < w_V + D0/tau_w and V <= V_T'
EQ_MEMBR_EVENT['w_crossing']['vars'] = ["V", "w"]
EQ_MEMBR_EVENT['w_crossing']['reset'] = "w=w_V - D0/tau_w"

IRHEOBASE = lambda g_L, V_T, E_L ,delta_T: g_L * (V_T - E_L - delta_T)


EQ_SYN_CHANNELS = '\n'.join(['{0}: 1'.format(name) for name in CHANNEL_NAMES])
EQ_STSP_VARS = '\n'.join(['{}: {}'.format(var, UNIT_MAIN_DICT[STSP_DECL_VARS[var]['unit']]) for var in STSP_DECL_VAR_NAMES])
EQ_SYN_PARAMS = '\n'.join(['{}: {}'.format(par, UNIT_MAIN_DICT[SYN_PARAMS[par]['unit']]) for par in SYN_PARAM_NAMES])
EQ_SYN_AUX = '\n'.join(['{}: {}'.format(par, UNIT_MAIN_DICT[SYN_AUX[par]['unit']]) for par in SYN_AUX_NAMES])

EQ_SYN_MODEL = '\n\n'.join([EQ_SYN_CHANNELS,EQ_STSP_VARS,
                            EQ_SYN_PARAMS, EQ_SYN_AUX])

EQ_EXT_SYN_MODEL = '\n\n'.join([EQ_SYN_CHANNELS, EQ_SYN_PARAMS,EQ_SYN_AUX])

EQ_SYN_PATHWAY = []
EQ_SYN_PATHWAY.append(dict(eq='failure = int(rand()<pfail)',order=0, delay=False))
EQ_SYN_PATHWAY.append(dict(eq="u_temp = U + u * (1 - U) * exp(-(t - last_spike_pre)/tau_fac)", order=1, delay=False))
EQ_SYN_PATHWAY.append(dict(eq="R_temp = 1 + (R - u * R - 1) * exp(- (t - last_spike_pre)/tau_rec)", order=2, delay=False))
EQ_SYN_PATHWAY.append(dict(eq='u = u_temp', order=3, delay=False))
EQ_SYN_PATHWAY.append(dict(eq='R = R_temp', order=4, delay=False))
EQ_SYN_PATHWAY.append(dict(eq='last_spike_pre = t', order=5, delay=False))
EQ_SYN_PATHWAY.append(dict(eq='a_syn = u * R', order=6, delay=False))

for name in CHANNEL_NAMES:    
    EQ_SYN_PATHWAY.append(dict(eq="g_{0}_on_post += {0} * gmax * a_syn * (1 - failure) * gsyn_amp".format(name), order=7, delay=True))
    EQ_SYN_PATHWAY.append(dict(eq="g_{0}_off_post += {0} * gmax * a_syn * (1 - failure) * gsyn_amp".format(name), order=7, delay=True))


EQ_EXT_SYN_PATHWAY = []
EQ_EXT_SYN_PATHWAY.append(dict(eq='failure = int(rand()<pfail)',order=0, delay=False))
for name in CHANNEL_NAMES:    
    EQ_EXT_SYN_PATHWAY.append(dict(eq="g_{0}_on_post += {0} * gmax  * (1 - failure) * gsyn_amp".format(name), order=0, delay=False))
    EQ_EXT_SYN_PATHWAY.append(dict(eq="g_{0}_off_post += {0} * gmax * (1 - failure) * gsyn_amp".format(name), order=0, delay=False))

EQ_VAR_UNITS = dict(V ='mV', w='pA', t='ms')
EQ_VAR_UNITS = EQ_VAR_UNITS | dict(I_tot='pA', I_syn='pA', I_inj='pA', I_DC='pA', I_AC='pA')
for name in CHANNEL_NAMES:
    EQ_VAR_UNITS['I_{0}'.format(name)]='pA'
for VAR_UNITS in [MEMB_PARAM_UNITS, STSP_DECL_VARS, SYN_PARAMS, SYN_AUX]:
    for var in VAR_UNITS:
        EQ_VAR_UNITS[var]=VAR_UNITS[var]['unit']



@time_report('Basics setup')
def basics_setup(Ncells_prompted, Nstripes, basics_scales=None, disp=True):
    
    # Copies of basic_scales targets are necessary when many simulations are carried out without
    # re-importing this module; otherwise, the modifications by basics_scales would remain
    # in the DataArrays between subsequent simulations
    PCON_COPY = PCON.copy()
    CELL_STD_COPY = CELL_STD.copy()
    SYN_GMAX_COPY = SYN_GMAX.copy()
    spiking_params_copy = spiking_params.copy()
    spiking_params_copy['gmax'] = dict(mean=SYN_GMAX_COPY, sigma=GMAX_SIGMA, min=GMAX_MIN, max=GMAX_MAX)

    
    SCALABLES = {'pCon': PCON_COPY, 'membr_param_std': CELL_STD_COPY, 'gmax_mean': SYN_GMAX_COPY}

    if basics_scales is not None:
    
        for param in basics_scales:              
            for targetsource_dict, scale in basics_scales[param]:
                SCALABLES[param].loc[targetsource_dict] = SCALABLES[param].loc[targetsource_dict].values * scale
                
    group_setup = GroupSetup(GROUP_KINDS, GROUPS_SETS)          
    stripes_setup = StripeSetup(Nstripes, INTERSTRIPE_SETS)
    connection_setup = ConnectionSetup(PCON_COPY, CLUSTER_FLAG)          
    cellsetup = StructureSetup(Ncells_prompted, PCELLS_PER_GROUP, group_setup, stripes_setup, connection_setup)
    
    membrane_setup = MembraneSetup(MEMB_PARAM_NAMES, MEMB_PARAM_UNITS, UNITS_DICT, UNIT_MAIN_DICT,
                                   CELL_MEAN, CELL_COVARIANCE, K_TRANSF, 
                                   CELL_STD_COPY, CELL_MINIMUM, CELL_MAXIMUM, TAU_MIN, TAU_MAX) 
    
    STSP_setup=STSPSetup(decl_vars=STSP_DECL_VARS,kinds=STPS_KINDS, sets=STSP_SETS, distr=STSP_SET_DISTR)
    channels_setup = ChannelSetup(CHANNEL_KINDS,PARAM_CHANNELS, GSYN_FACTOR_CHANNELS)      
    
    spike_params_data = ParamSetup()
    for par in SYN_PARAM_NAMES:
        spike_params_data[par] = spiking_params_copy[par]
        
    
    
    # spike_params_setup = SpikeParamSetup(SYN_PARAMS, SYN_PFAIL, SYN_GMAX, GMAX_SIGMA, GMAX_MIN, GMAX_MAX)
    spike_params_setup = SpikeParamSetup(SYN_PARAMS, spike_params_data)
    delay_setup = DelaySetup(SYN_DELAY, DELAY_SIGMA, DELAY_MIN, DELAY_MAX)
    synapses_setup = SynapseSetup(SYN_KINDS, spike_params_setup, delay_setup, STSP_setup, channels_setup)
    equations_setup = EquationsSetup(EQ_MEMBR_MODEL,EQ_MEMBR_THRESHOLD, EQ_MEMBR_RESET, EQ_MEMBR_EVENT, IRHEOBASE,
                                     EQ_SYN_MODEL, EQ_SYN_PATHWAY, EQ_EXT_SYN_MODEL, EQ_EXT_SYN_PATHWAY,
                                     EQ_VAR_UNITS,)
    
    

    if cellsetup.Ncells != Ncells_prompted and disp:
        print('REPORT: The number of neurons was adjusted from {} to {} to preserve group distribution.\n'.format(Ncells_prompted, cellsetup.Ncells))
      
    return BasicsSetup(cellsetup, membrane_setup, synapses_setup, equations_setup, SCALABLES, basics_scales)


def setup_view(setup):
    try:
        setup.view()
    except AttributeError:
        print('Input is {}, not BaseClass.'.format(setup.__class__.__name__))
        print(setup)
def setup_tree(setup):
    try:
        setup.tree()
    except AttributeError:
        print('Input is {}, not BaseClass.'.format(setup.__class__.__name__))
        print(setup)




@dataclass
class BaseClass:
    pass

    def view(self, item=None):
        self.tree(0, item)
    
    def tree(self, max_depth=None, item=None):
        if item is None:
            inst=self
        else:
            inst=self[item]
      
        
        if len(inst.values())==0:
            print('{}'.format(inst.__class__.__name__))
        else:
            print('{} instance containing:'.format(inst.__class__.__name__))       
            for key in  inst.keys():
                if getattr(inst[key], 'subtree', None) is not None and (max_depth is None or max_depth>0):
                    post = ' instance containing:'
                else:
                    post = ''
                print(' '*4+'|- {}: {}'.format(key, inst[key].__class__.__name__)+post)
                # print(max_depth)
                if getattr(inst[key], 'subtree', None) is not None:
                    if max_depth is None or max_depth>0:
                        inst[key].subtree(1, max_depth)    
    
    def subtree(self, depth, max_depth):
        tab = '    |' *depth
        # post = ''+' instance containing:'
        # print(tab + tab1 + '{} instance containing:'.format(self.__class__.__name__))       
        for key in  self.keys():
            if getattr(self[key], 'subtree', None) is not None and (max_depth is None or max_depth>depth):
                post = ' instance containing:'
            else:
                post = ''
            print(tab+ ' '*4+'|- {}: {}'.format(key, self[key].__class__.__name__)+post)
            if getattr(self[key], 'subtree', None) is not None:
                if max_depth is None or max_depth>depth:
                    self[key].subtree(depth+1, max_depth)

    def keys(self):
        return list(self.__dict__.keys())
    
    
    def values(self):
        return list(self.__dict__.values())
    
    def __getitem__(self, item):
        return getattr(self, item)
    
    def __setitem__(self, item, value):
        setattr(self, item, value)
    
    
@dataclass
class ChannelSetup(BaseClass):
    kinds: list[str]
    params: xr.DataArray
    gsyn_factor: xr.DataArray
    names: list = field(default_factory=list)
    kinds_to_names: dict = field(default_factory=dict)
    
    def __post_init__(self):
        self.names = list(self.kinds.keys())
        self.kinds_to_names = {'exc':[], 'inh':[]}
        for name in self.names:
            self.kinds_to_names[self.kinds[name]].append(name)


@dataclass
class SynapseSetup(BaseClass):
    kinds: any
    spiking: any
    delay: any
    STSP: any
    channels: any

@dataclass
class SpikeParamSetup(BaseClass):
    names: any
    params: any
    
@dataclass
class DelaySetup(BaseClass):
    delay_mean: any
    delay_sigma: any
    delay_min: any
    delay_max: any
    

@dataclass
class StripeSetup(BaseClass):
    N: int
    inter: any


@dataclass
class ConnectionSetup(BaseClass):
    pCon: any
    cluster: any
    

@dataclass
class STSPSetup(BaseClass):
    decl_vars: any
    kinds: dict
    sets: any
    distr: dict

@dataclass
class MembraneSetup(BaseClass):  
    names: any
    name_units: any
    unit_dict: any
    unit_main_dict: any
    mean: float
    covariance: float
    k: float
    std: float
    min: float
    max:float
    tau_m_min: any
    tau_m_max: any

 
@dataclass
class GroupSetup(BaseClass):   
    kinds: any
    sets:any
    names: list[str] = field(default_factory=list)
    N:int = 0
    
    def __post_init__(self):
        self.names = list(self.kinds.keys())
        self.N = len(self.names)
        self.idcs = {}
        for name_idc in range(len(self.names)):
            self.idcs[self.names[name_idc]] = name_idc


@dataclass
class StructureSetup(BaseClass):
         
    Ncells_prompt: int
    Pcells_per_group: list or np.array
    groups: GroupSetup
    stripes: dict
    conn: any
    Ncells: int = 0
    Ncells_total: int = 0
    Ncells_per_group: xr.DataArray = xr.DataArray([])
        
    def __post_init__(self):
        self.Ncells_per_group = xr.DataArray(np.ceil((self.Ncells_prompt*self.Pcells_per_group)/100).astype(int), coords=[self.groups.names], dims='group')
        self.Ncells = int(sum(self.Ncells_per_group))
        self.Ncells_total = self.stripes.N*self.Ncells
        #self.groups = GroupSetup(names, sets, types)
        
    
@dataclass
class ParamSetup(BaseClass):
    pass

@dataclass
class EquationsSetup(BaseClass):
    membr_model: str
    membr_threshold: str
    membr_reset: str
    membr_events: dict
    rheobase: any
    syn_model: str
    syn_pathway: list
    ext_syn_model: str
    ext_syn_pathway: list
    var_units: dict
        
@dataclass
class BasicsSetup(BaseClass):
    struct: any
    membr: any
    syn: any
    equations: any
    scalables: any
    scales: any



if __name__ == '__main__':
    
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
    N1=control_param['Neurons per stripe']
    Nstripes=control_param['Stripes']
    gmax_scale_prompted = [1, 1, 1, 1]
    pCon_scale = [1, 1, 1, 1]
    param_std_scale = np.ones((9, 14))

    # scales = gmax_scale_prompted, pCon_scale, param_std_scale
    gmax_scales = [(dict(target=GROUPS_SETS['PC'], source=GROUPS_SETS['PC']), 2)]
    basics = basics_setup(N1, Nstripes)
    
    

    
