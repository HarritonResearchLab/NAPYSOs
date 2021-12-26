# Import(s)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Definitions

def cody_Q(mjd, mag, magerr, timescale, sig_factor):
    # Import(s)
    from astropy.convolution import Box1DKernel, convolve
    
    # Convert to arrays
    mjd = np.array(mjd)
    mag = np.array(mag)
    magerr = np.array(magerr)
        
    # Calculate sig
    sig = sig_factor*np.mean(magerr)

    # Create the residual curve
    phase = mjd % timescale
    mag = mag[np.argsort(phase)]

    # We use three periods and extract the middle to prevent edge effects
    three_periods = np.concatenate((mag, mag, mag))
    boxcar = Box1DKernel(len(mag) // 4)
    smooth_mag = convolve(three_periods, boxcar)
    smooth_mag = smooth_mag[np.size(mag):2*np.size(mag)]

    resid_mag = mag - smooth_mag

    q = float((np.nanstd(resid_mag) ** 2 - sig ** 2)/(np.nanstd(mag) ** 2 - sig ** 2))

    return q, resid_mag, smooth_mag 

data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
id = 'FHK_26'
per = 13.704636645102322
    
data_df = pd.read_csv(data_temp.replace('+++', id))
mjd = np.array(data_df['mjd'])
mag = np.array(data_df['mag'])
magerr = np.array(data_df['magerr'])

smoothed_mag = cody_Q(mjd, mag, magerr, per, 1.25)[2]
# Import(s)
import os 

# Action
plt.rcParams['font.family']='serif'

plt.tick_params(axis='both',which='both',direction='in')
plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)

# Our period 
phased_dates = (np.mod(mjd, per))/per 
phased_dates_cycle_2 = phased_dates + 1

plt.scatter(phased_dates, mag, color='#408ee0', s=3)
plt.scatter(phased_dates_cycle_2, mag, color='#408ee0', s=3)
plt.scatter(phased_dates, smoothed_mag, color='red')
plt.scatter(phased_dates_cycle_2, smoothed_mag, color='red')
plt.xlabel('Our Folded LC (P='+str(round(per, 2))+')')
plt.ylabel('r (mag)')

plt.show()