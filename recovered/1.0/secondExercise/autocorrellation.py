import matplotlib.pyplot as plt 
from astropy.io import ascii
import numpy as np
import statsmodels.api as sm
from personalastropy.caltech.Caltech import sort_data

#get data
green_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_g.tbl'
red_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_r.tbl'
g = ascii.read(green_path,format='ipac',delimiter='|')
r = ascii.read(red_path,format='ipac',delimiter='|')

green_dates = np.array(g['hjd'])
green_mags = np.array(g['mag'])
red_dates = np.array(r['hjd'])
red_mags = np.array(r['mag'])

#sort data
sgd = sort_data(green_dates,green_mags)[0]
sgm = sort_data(green_dates,green_mags)[1]
srd = sort_data(red_dates,red_mags)[0]
srm = sort_data(red_dates,red_mags)[1]

import math
x = np.linspace(0,18,40)
y = np.sin(x)
x2 = x-1
y2 = np.sin(x2)

plt.scatter(y,y2)


plt.show()

