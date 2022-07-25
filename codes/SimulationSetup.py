import numpy as np
import os

constant_current = [[250, 200, 200, 200, 200, 200, 200, 250, 200, 200, 200, 200, 200, 200],
                    ]
fluctuating_current= []
PoissonStimuli = list()
RegularStimuli = list()
NeuronMonitor = list()
SynapsesMonitor = list()
analysis_params = dict()
analysis_params['ISI analysis'] = dict()
analysis_params['ISI correlation1'] = dict()
analysis_params['ISI correlation2'] = dict()
analysis_params['ISI Pearson'] = dict()
analysis_params['V analysis'] = dict()
analysis_params['V correlation'] = dict()
analysis_params['Frequence'] = dict()
analysis_params['Populational rate'] = dict()
analysis_params['Rate stratification'] = dict()


def set_dir():
    
    k = 0
    simulation_dir = 'Simulation_{}'.format(k)
    while os.path.isdir(simulation_dir):
        k+=1
        simulation_dir = 'Simulation_{}'.format(k)
        
    os.mkdir(simulation_dir)
    
    return simulation_dir


def save_group_distr(group_distr, sim_dir):
    
    with open('{}/group_distr.txt'.format(sim_dir), 'a') as f:
        for ii in range(len(group_distr[0])):
            for i in range(len(group_distr)):
                print(*group_distr[i][ii], sep=',', end='',file=f)
                if i != len(group_distr)-1:
                    print('-', end='', file=f)
            if ii != len(group_distr[0])-1:
                print(';', end='', file=f)
     
                
def set_group_distr(sim_dir):
    
    with open('Simulation_{}/group_distr.txt'.format(sim_dir), 'r') as f:
        tex = f.read()

    stripes = tex.split(';')
    Nstripes = len(stripes)
    Ngroups = len(stripes[0].split('-'))

    distr = [[[] for j in range(Nstripes)] for i in range(Ngroups)]

    for j in range(Nstripes):        
        stripe = tex.split(';')[j]
        
        for k in range(Ngroups):    
            group = stripe.split('-')[k]
            if group != '':
                distr[k][j].extend(np.asarray(group.split(',')).astype(int))
            
    return distr


def save_setup(NeuPar, V0, STypPar, SynPar,SPMtx, group_distr, sim_dir):
    
    np.savetxt('{}/NeuPar.txt'.format(sim_dir), NeuPar, delimiter=',')
    np.savetxt('{}/V0.txt'.format(sim_dir), V0, delimiter=',')
    np.savetxt('{}/STypPar.txt'.format(sim_dir), STypPar, delimiter=',')
    np.savetxt('{}/SynPar.txt'.format(sim_dir), SynPar, delimiter=',')
    np.savetxt('{}/SPMtx0.txt'.format(sim_dir), SPMtx[:, :, 0], delimiter=',')
    np.savetxt('{}/SPMtx1.txt'.format(sim_dir), SPMtx[:, :, 1], delimiter=',')
    
    save_group_distr(group_distr, sim_dir)
    
    
   
    
def set_setup(sim_dir):
    NeuPar = np.genfromtxt('Simulation_{}/NeuPar.txt'.format(sim_dir), delimiter=',')
    V0 = np.genfromtxt('Simulation_{}/V0.txt'.format(sim_dir), delimiter=',')
    STypPar = np.genfromtxt('Simulation_{}/STypPar.txt'.format(sim_dir), delimiter=',')
    SynPar = np.genfromtxt('Simulation_{}/SynPar.txt'.format(sim_dir), delimiter=',')
    SPMtx0 = np.genfromtxt('Simulation_{}/SPMtx0.txt'.format(sim_dir), delimiter=',')
    SPMtx1 = np.genfromtxt('Simulation_{}/SPMtx1.txt'.format(sim_dir), delimiter=',')
    
    SPMtx = np.zeros((*SPMtx0.shape, 2))
    SPMtx[:, :, 0] = SPMtx0
    SPMtx[:, :, 1] = SPMtx1
    
    group_distr = set_group_distr(sim_dir)
    
    return NeuPar, V0, STypPar, SynPar, SPMtx, group_distr


def save_info(control_param, clustering, PoissonStimuli,
              RegularStimuli, NeuronMonitor, SynapsesMonitor, analysis_params,
              scales1, scales2, constant_current, fluctuating_current, sim_dir, file='SIMULATION_INFO'):
    
    with open('{}/{}_{}.txt'.format(sim_dir, file, sim_dir), 'a') as f:
        
        print('#'*20, 'SIMULATION INFO', '#'*20, file=f)
        
        control_param_str = '{\n'
        
        for key, value in control_param.items():
            control_param_str += "'{}': {},\n".format(key, value)
        control_param_str += '}'
                
        
        clustering_str = 'clustering: {}'.format(clustering)
        
        analysis_params_str = '{\n'
        
        for key, value in analysis_params.items():
            analysis_params_str += "'{}': {},\n".format(key, value)
        analysis_params_str += '}'

        PoissonStimuli_str = '[\n'
        
        for item in PoissonStimuli:
            PoissonStimuli_str += "{},\n".format(item)
        PoissonStimuli_str += ']'
        
        RegularStimuli_str = '[\n'
        
        for item in RegularStimuli:
            RegularStimuli_str += "{},\n".format(item)
        RegularStimuli_str += ']'
        
        NeuronMonitor_str = '[\n'
        
        for item in NeuronMonitor:
            NeuronMonitor_str += "{},\n".format(item)
        NeuronMonitor_str += ']'
        
        SynapsesMonitor_str = '[\n'
        
        for item in SynapsesMonitor:
            SynapsesMonitor_str += "{},\n".format(item)
        SynapsesMonitor_str += ']'
        
        constant_current_str = '[\n'
        
        for item in constant_current:
            constant_current_str += "{},".format(item)
        constant_current_str += '\n]'
        
        print('control_param\n', control_param_str, end='\n'+'-'*50+'\n', file=f)
        print('clustering\n', clustering_str, end='\n'+'-'*50+'\n', file=f)
        print('PoissonStimuli\n', PoissonStimuli_str, end='\n'+'-'*50+'\n', file=f)
        print('RegularStimuli\n', RegularStimuli_str, end='\n'+'-'*50+'\n', file=f)
       
        print('NeuronMonitor\n', NeuronMonitor_str, end='\n'+'-'*50+'\n', file=f)
        print('SynapsesMonitor\n', SynapsesMonitor_str, end='\n'+'-'*50+'\n', file=f)
        print('analysis_params\n', analysis_params_str, end='\n'+'-'*50+'\n', file=f)
        print('constant_current\n', constant_current_str, end='\n'+'-'*50+'\n', file=f)
        print('scales1\n', scales1, end='\n'+'-'*50+'\n', file=f)
        print('scales2\n', scales2, end='\n'+'-'*50+'\n', file=f)
        print('fluctuating_current\n', fluctuating_current, end='\n'+'-'*50+'\n', file=f)