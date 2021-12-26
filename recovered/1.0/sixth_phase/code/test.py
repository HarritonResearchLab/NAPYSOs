#Import(s)
import numpy as np
from astropy.io import ascii
import math
import matplotlib.pyplot as plt

#Action

#Some defs
r = ascii.read('/home/thaddaeus/FMU/HRL/LAH/fourth_phase/data/extended_lightcurves/lc_20:05:06.4+36:29:13.3_r.tbl',format='ipac',delimiter='|')
r_dates = np.array(r['mjd'])
r_mags = list(r['mag'])
r_magerrs = list(r['magerr'])
tolerance=4
paired_arrays = [r_mags,r_magerrs]

#quick plot
fig, axs = plt.subplots(2,1,sharex=True)
axs[0].scatter(r_dates,r_mags,color='indianred',marker='.',s=0.5)
axs[0].invert_yaxis()
axs[0].set_xlabel("Time (MJD)")
axs[0].set_ylabel("Mag (r)")

#Actual action
min_date = math.floor(np.min(r_dates))
max_date = math.ceil(np.max(r_dates))
date_bins = list(range(min_date,max_date,1))
binned_dates = np.digitize(x=r_dates,bins=date_bins)
bin_counts = np.bincount(binned_dates)

nz_bin_counts = bin_counts[np.nonzero(bin_counts)]

mean_counts = np.mean(nz_bin_counts)
std_counts = np.std(nz_bin_counts,ddof=1)

cutoff = mean_counts+tolerance*(std_counts)

intervals_to_remove = []
bin_counts = list(bin_counts)

for count in bin_counts:
    if count > cutoff:
        bin_index = bin_counts.index(count)
        lower_edge = date_bins[bin_index-1]
        upper_edge = date_bins[bin_index]
        edge_string = str(lower_edge)+':'+str(upper_edge)
        intervals_to_remove.append(edge_string)

indices_to_remove = []

r_dates = list(r_dates)

for interval in intervals_to_remove:
    interval_list = interval.split(':')
    lower = float(interval_list[0])
    upper = float(interval_list[1])
    for date in r_dates:
        if date > lower and date < upper:
           date_index = r_dates.index(date)
           indices_to_remove.append(date_index)

for i in sorted(indices_to_remove, reverse=True):
    del r_dates[i]
    for paired_array in paired_arrays:
        del paired_array[i]





axs[1].scatter(r_dates,r_mags,color='indianred',marker='.',s=0.5)
axs[1].invert_yaxis()
axs[1].set_xlabel("Time (MJD)")
axs[1].set_ylabel("Mag (r)")

plt.show()




