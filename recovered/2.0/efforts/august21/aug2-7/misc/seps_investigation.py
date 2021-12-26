import numpy as np

sep_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/zeta/prepping_stuff/seps.txt'
seps = np.loadtxt(fname=sep_file)

def return_percentile(seps_arr, sep):
    count = len(seps[np.where(seps_arr<=sep)])
    percentile = count/len(seps_arr)
    return percentile
    
    
print('0.25: ', return_percentile(seps, 0.25), 
'1.25: ',return_percentile(seps, 1.25))