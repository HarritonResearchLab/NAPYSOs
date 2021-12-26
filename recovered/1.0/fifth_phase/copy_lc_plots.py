#Import(s)
import numpy as np
import os 
import shutil

from scipy.stats import kurtosis



coords = []
with open('/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/coordinates.txt','r') as f:
    for line in f:
        line = line.replace('\n','')
        coords.append(line)

for coord in coords:
    old_loc = '/home/thaddaeus/FMU/HRL/LAH/third_phase/images/lightcurve_plots/+++++.png'.replace('+++++',coord)
    if os.path.exists(old_loc)==True:
        new_loc = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/for_email/plots/lc_plots/+++++.png'.replace('+++++',coord)
        shutil.copy(old_loc,new_loc)
