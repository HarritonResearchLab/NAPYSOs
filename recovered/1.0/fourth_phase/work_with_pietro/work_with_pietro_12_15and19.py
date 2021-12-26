from astropy.io import ascii 
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
import pandas as pd

#####CHANGE THESE\/\/\/\/######

obsids_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/key.csv'
r_path_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++_r.tbl'
g_path_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++_g.tbl'
#############

coordinates = list(pd.read_csv(obsids_file)['ID'])



r_median_mags = []
r_median_errors = []
g_median_mags = []
g_median_errors = []

for item in coordinates:
    r_path = r_path_temp.replace('+++',item)
    r = ascii.read(r_path,format='ipac',delimiter='|')
    red_dates = np.array(r['hjd'])
    red_mags = np.array(r['mag'])
    red_mag_errors = np.array(r['magerr'])

    g_path = g_path_temp.replace('+++',item)
    g = ascii.read(g_path,format='ipac',delimiter='|')
    green_dates = np.array(g['hjd'])
    green_mags = np.array(g['mag'])
    green_mag_errors = np.array(g['magerr'])

    if len(red_dates) > 3 and len(red_mags) > 3:
        median_mag = np.median(red_mags)
        r_median_mags.append(median_mag)
        median_error = np.median(red_mag_errors)
        r_median_errors.append(median_error)
    
    if len(green_dates) > 3 and len(green_mags) > 3:
        median_mag = np.median(green_mags)
        g_median_mags.append(median_mag)
        median_error = np.median(green_mag_errors)
        g_median_errors.append(median_error)
        
zipped_list = list(zip(g_median_mags,g_median_errors,r_median_mags,r_median_errors))
df = pd.DataFrame(data=zipped_list,columns=['g_median_mag', 'g_median_magerr', 'r_median_mag', 'r_median_magerr'])
df.to_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/alpha/for_pietro/data.csv',index=False)



def exponential(x,a,b,c): 
    return((a*np.exp(b*x))+c)

pars, cov = curve_fit(f=exponential, xdata=g_median_mags, ydata=g_median_errors)

print("Formula: a*b^(x)+c")
print('Pars: a,b,c')
print(pars)

best_fit_ys = []

for item in list(g_median_mags):
    fitted = (pars[0]*np.exp(pars[1]*item))+pars[2]
    best_fit_ys.append(fitted)
best_fit_ys = np.array(best_fit_ys)

plt.scatter(g_median_mags,g_median_errors)
plt.scatter(g_median_mags,best_fit_ys,label='Best Exponential Fit',color='orange',s=1)

plt.xlabel("Median Mag")
plt.ylabel('Median Magerr')
plt.legend()
plt.title('Green')
plt.show()
