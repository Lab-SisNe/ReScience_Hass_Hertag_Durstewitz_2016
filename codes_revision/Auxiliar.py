import os
import time as tm

def set_simulation_dir(name=None):
    if name is None:
        num = 0
        name = 'Simulation_{}'
        while os.path.isdir(name.format(num)):
            num+=1
        name=name.format(num)
    elif os.path.exists(name):
        raise NameError('Path already exists.')
    
    os.mkdir(name)
    print('REPORT: simulation_dir created:', name, end='\n\n')
    
    return name

def time_report(text=None):
    def wrap(f):    
        def wrapped_f(*args, **kwargs):
            start_time = tm.time()
            print_out1 = '|   TIME REPORT:   |'
            print_out2 =  ' Starting {} ...'.format(text if text is not None else f.__name__)
            print('-'*len(print_out1) + ' '*len(print_out2))
            print(print_out1+print_out2)
            print('-'*len(print_out1) + ' '*len(print_out2), end='\n\n')
            ret = f(*args, **kwargs)
            elapsed = tm.time()-start_time
            if elapsed>=3600:
                print_out = '||   TIME REPORT: {} completed in {:.0f} h, {:.0f} min and {:.2f} s.   ||'.format(text if text is not None else f.__name__, elapsed//3600, elapsed//60, elapsed%60)
                print('-'*len(print_out))
                print(print_out)
            elif elapsed>=60:
                print_out = '||   TIME REPORT: {} completed in {:.0f} min and {:.2f} s.   ||'.format(text if text is not None else f.__name__, elapsed//60, elapsed%60)
                print('-'*len(print_out))
                print(print_out)
            else:
                print_out = '||   TIME REPORT: {} completed in {:.2f} s.   ||'.format(text if text is not None else f.__name__, elapsed)
                print('-'*len(print_out))
                print(print_out)
            print('-'*len(print_out), end='\n\n')
            return ret
        return wrapped_f
    return wrap

