'''
import re 
import os, shutil
import requests
import gzip
import astropy 
import astroquery
import sys
import statistics as stat
import os.path
import numpy as np
import pandas as pd
from astroquery.mearth import Mearth
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
'''

def removeIntervals(x,y,intervals):
    #Import(s)
    import numpy as np 

    #Action
    dates = list(x)
    mags = list(y)
    
    for item in intervals:
        lower_bound = float(item.split(':')[0])
        upper_bound = float(item.split(':')[1])
        for elem in dates:
            if elem > lower_bound:
                if elem < upper_bound:
                    elem_index = dates.index(elem)
                    dates.remove(elem)
                    mags.remove(mags[elem_index])


    
    return np.array(dates),np.array(mags)


import numpy as np
from astropy.io import ascii
from personalastropy.caltech.Caltech import sort_data
import matplotlib.pyplot as plt 

red_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_r.tbl'
r = ascii.read(red_path,format='ipac',delimiter='|')
red_dates = np.array(r['hjd'])
red_mags = np.array(r['mag'])
srd = sort_data(unsorteddates=red_dates,unsortedmags=red_mags)[0]
srm = sort_data(unsorteddates=red_dates,unsortedmags=red_mags)[1]

new_srd = removeIntervals(srd,srm,intervals=['2458437:2458438','2458455:2458456'])[0]
new_srm = removeIntervals(srd,srm,intervals=['2458437:2458438','2458455:2458456'])[1]


plt.scatter(new_srd,new_srm,c='red',s=2)
plt.xlabel('HJD')
plt.ylabel('Mag')
plt.gca().invert_yaxis()
plt.show()
    

'''
num_arrays = returnGoodRegions(x=sgd,y=sgm,max_sep=20,min_card=3)[0]

import matplotlib.pyplot as plt 

x_arrays = returnGoodRegions(x=sgd,y=sgm,max_sep=20,min_card=3)[1]
y_arrays = returnGoodRegions(x=sgd,y=sgm,max_sep=20,min_card=3)[2]
i = 0
colors = ['red','blue','green','yellow','purple','black','brown','pink','orange']
for elem in x_arrays:
    plt.scatter(elem,y_arrays[i],s=4,c=colors[i])
    i += 1
plt.gca().invert_yaxis()
plt.xlabel('HJD')
plt.ylabel('Mag')
plt.show()
'''

    
