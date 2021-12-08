#Import(s)
import numpy as np
import os 
import shutil
from astropy.io import ascii
from scipy.stats import kurtosis



coords = []
with open('/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/coordinates.txt','r') as f:
    for line in f:
        line = line.replace('\n','')
        coords.append(line)

for coord in coords:
    r_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/lightcurves/lc_+++++_r.tbl'.replace('+++++',coord)

    r = ascii.read(r_file,format='ipac',delimiter='|')

    red_mags=np.array(r['mag'])
    red_dates=np.array(r['hjd'])
    red_magerrs=np.array(r['magerr'])

    red_std = np.std(red_mags,ddof=1)

    red_kurt = kurtosis(red_mags)
    
    if red_std > 0.6:
        print(coord)
    
    if red_kurt > 30:
        print(coord)
    
