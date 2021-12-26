from astropy.io import ascii 
import numpy as np
import matplotlib.pyplot as plt
import os

coordinates = []
with open('/home/thaddaeus/FMU/HRL/LAH/third_phase/data/obsids.txt','r') as f:
    for line in f:
        line = line.replace('\n','')
        if os.path.exists('/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_'+line+'_r.tbl'):
            coordinates.append(line)

median_mags = []
median_errors = []
stds = []

for item in coordinates:
    path = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_'+item+'_r.tbl' 
    r = ascii.read(path,format='ipac',delimiter='|')
    red_dates = np.array(r['hjd'])
    red_mags = np.array(r['mag'])
    red_mag_errors = np.array(r['magerr'])

    if len(red_dates) > 3 and len(red_mags) > 3:
        median_mag = np.median(red_mags)
        median_mags.append(median_mag)
        median_error = np.median(red_mag_errors)
        median_errors.append(median_error)
        stdev = np.std(red_mags,ddof=1)
        stds.append(stdev)

plt.scatter(median_mags,median_errors)
plt.xlabel("Median Mag")
plt.ylabel('Median Magerr')
plt.show()