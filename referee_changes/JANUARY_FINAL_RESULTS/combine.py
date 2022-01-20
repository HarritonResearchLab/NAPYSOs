from asyncio import format_helpers
import pandas as pd
import numpy as np


### CALC Q FROM CODE DIR ####

def Q(dates, mags, magerrs, timescale, sig_factor): 
    '''
    Parameters: 
    dates: array of time values associated with observations
    mags: array of magnitude values associated with observations
    magerrs: array of errors associated with magnitudes in mags array
    timescale: float representative of an object's timescale of variability
    sig_factor: multiplied against mean magerr, should be set to 1 initially
    
    Returns: 
    q: quasi-periodicity value
    resid_mags: array of residual mags
    '''
    
    # Import(s)
    import numpy as np
    from astropy.convolution import Box1DKernel, convolve
    
    # Convert to arrays
    mjd = np.array(dates)
    mag = np.array(mags)
    magerr = np.array(magerrs)
        
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

    return q, resid_mag


#### ACTUAL WORK ####


september_df = pd.read_csv(r'C:\Users\Research\Documents\GitHub\NAPYSOs\recovered\2.0\efforts\september21\SEPTEMBER_FINAL_RESULTS\merged_results.csv')
jan_df = pd.read_csv(r'C:\Users\Research\Documents\GitHub\NAPYSOs\referee_changes\lah_vetting_results\ref_files\lah_changes_incorporated.csv')

merged_df = september_df.merge(jan_df, on='ID')
merged_df = merged_df.drop(labels=['Q', 'M', 'primary_class', 'secondary_class', 'RA_jan', 'DEC_jan', 'PER'], axis='columns')
merged_df = merged_df.rename(columns={'Q_jan':'Q', 'M_jan':'M', 'primary_class_jan':'primary_class', 'secondary_class_jan':'secondary_class', 'PER_jan':'PER'})

qs = np.array(merged_df['Q'])
pers = np.array(merged_df['PER'])
changed_pers = np.array(merged_df['per_changed'])
ids = np.array(merged_df['ID'])

for index, i in enumerate(changed_pers): 
    if type(i)!=float: 
        if i==True: 
            lc_file = r"./recovered/2.0/data/AUGUST_5th/light_curves/" 
            lc_file += ids[index] + '_r.csv'
            lc_df = pd.read_csv(lc_file)
            mags = np.array(lc_df['mag'])
            dates = np.array(lc_df['mjd'])
            magerrs = np.array(lc_df['magerr'])

            new_q = Q(dates, mags, magerrs, pers[index], 1.25)[0]

            if 0 > new_q or new_q > 1: 
                print(new_q, ids[index])
            
            qs[index] = new_q

merged_df['Q'] = qs

merged_df = merged_df.drop(columns=['odr slope','odr slope error','intercept','intercept error','angle','angle error','pythag'])

odr_df = pd.read_csv('https://raw.githubusercontent.com/HarritonResearchLab/NAPYSOs/main/referee_changes/odr-redux/odr_results.csv')

final_merged = merged_df.merge(odr_df, left_on='ID', right_on='star Identifier')

final_merged = final_merged.drop(columns=['per_changed', 'star Identifier'])

print(final_merged)


import matplotlib.pyplot as plt

plt.hist(final_merged['Q'])
plt.show()

final_merged.to_csv("./referee_changes/JANUARY_FINAL_RESULTS/JANUARY_FINAL_RESULTS.csv", index=False)