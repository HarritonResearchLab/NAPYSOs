#Import(s)
from personalastropy.ysospy.variability import lombScargle
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

#Action

red_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/FHK_286_r.tbl'
r = ascii.read(red_file,format='ipac',delimiter='|') 
o_red_dates = list(r['hjd'])
o_red_mags = list(r['mag'])


red_dates = []
for item in o_red_dates:
    item = float(item)
    red_dates.append(item)

red_mags = []
for item in o_red_mags:
    item = float(item)
    red_mags.append(item)


n_r_dates = []
n_r_mags = []
for date in red_dates:
    if date < 2458448 or date > 2458456:
        n_r_dates.append(date)
        date_index = red_dates.index(date)
        n_r_mags.append(red_mags[date_index])

#Plot
plt.scatter(n_r_dates,n_r_mags,color='indianred',marker='.')
plt.xlabel('Time (MJD)')
plt.ylabel('Mag')
plt.gca().invert_yaxis()
plt.show()


lombScargle(id='FHK_176',x=n_r_dates,y=n_r_mags,min_period=2,max_period=400,out_type='FHK_286_ls.png')