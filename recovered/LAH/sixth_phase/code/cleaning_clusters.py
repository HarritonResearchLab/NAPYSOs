#Import(s)
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt

#Action

r = ascii.read('/home/thaddaeus/FMU/HRL/LAH/fourth_phase/data/extended_lightcurves/lc_20:05:06.4+36:29:13.3_r.tbl',format='ipac',delimiter='|')
r_dates = list(r['mjd'])
r_mags = list(r['mag'])
r_magerrs = list(r['magerr'])

def clean_clusters(dates,paired_lists,tolerance):
    
    #Import(s)
    import numpy as np
    import math

    #Action

    min_date = math.floor(np.min(dates))
    max_date = math.ceil(np.max(dates))
    date_bins = list(range(min_date,max_date,1))
    binned_dates = np.digitize(x=dates,bins=date_bins)
    bin_counts = np.bincount(binned_dates)

    nz_bin_counts = bin_counts[np.nonzero(bin_counts)] #Get rid of bin_counts = 0 for calculations

    mean_counts = np.mean(nz_bin_counts)
    std_counts = np.std(nz_bin_counts,ddof=1)

    cutoff = mean_counts+float(tolerance)*(std_counts)

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

    for interval in intervals_to_remove:
        interval_list = interval.split(':')
        lower = float(interval_list[0])
        upper = float(interval_list[1])
        for date in dates:
            if date > lower and date < upper:
                date_index = dates.index(date)
                indices_to_remove.append(date_index)

    for i in sorted(indices_to_remove, reverse=True):
        del dates[i]
        for paired_list in paired_lists:
            del paired_list[i]

    return dates, paired_lists

print(clean_clusters(dates=r_dates,paired_lists=[r_mags,r_magerrs],tolerance=4)[1][1])



'''
r_dates =cleaner[0]
r_mags = cleaner[1][0]
r_magerrs = cleaner[1][1]
'''

def plot_it():
    plt.style.use('ggplot')
    plt.errorbar(r_dates,r_mags,yerr=r_magerrs,lw=0,elinewidth=0.5,marker='.',ms=1)
    plt.gca().invert_yaxis()
    plt.xlabel('Time (MJD)',fontsize=7)
    plt.ylabel('Mag',fontsize=7)

