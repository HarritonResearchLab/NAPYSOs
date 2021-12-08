#Import(s)
import numpy as np
from personalastropy.ysospy.variability import sokolovskyNu
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn as sns
#Action

coords = []
with open('/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/coordinates.txt','r') as f:
    for line in f:
        line = line.replace('\n','')
        coords.append(line)


ols_68_rsq = []
ols_80_rsq = []

ols_68_slopes = []
ols_80_slopes = []

ols_68_nu = []
ols_80_nu = []

ols_68_std = []
ols_80_std = []

ols_i_68_rsq = []
ols_i_80_rsq = []

ols_i_68_nu = []
ols_i_80_nu = []

ols_i_68_std = []
ols_i_80_std = []

ols_i_68_slopes = []
ols_i_80_slopes = []

for coord in coords:
    #std and nu
    
    r_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/lightcurves/lc_+++++_r.tbl'.replace('+++++',coord)

    r = ascii.read(r_file,format='ipac',delimiter='|')

    red_mags=np.array(r['mag'])
    red_dates=np.array(r['hjd'])
    red_magerrs=np.array(r['magerr'])

    red_std = np.std(red_mags,ddof=1)
    red_nu = sokolovskyNu(red_mags,red_magerrs)

    with open('/home/thaddaeus/FMU/HRL/LAH/fifth_phase/work/cmd_regression_results.csv','r') as f:
        for line in f:
            if coord in line:
                line_list = line.split(',')
                ols_slope = float(line_list[1])
                ols_rsq = float(line_list[2])
                ols_i_slope = float(line_list[3])
                ols_i_rsq = float(line_list[4])
                
                if ols_rsq > 0.68:
                    ols_68_rsq.append(ols_rsq)
                    ols_68_slopes.append(ols_slope)
                    ols_68_nu.append(red_nu)
                    ols_68_std.append(red_std)

                if ols_rsq > 0.8:
                    ols_80_rsq.append(ols_rsq)
                    ols_80_slopes.append(ols_slope)
                    ols_80_nu.append(red_nu)
                    ols_80_std.append(red_std)
                
                if ols_i_rsq > 0.68:
                    ols_i_68_rsq.append(ols_i_rsq)
                    ols_i_68_slopes.append(ols_i_slope)
                    ols_i_68_nu.append(red_nu)
                    ols_i_68_std.append(red_std)
                if ols_i_rsq > 0.8:
                    ols_i_80_rsq.append(ols_i_rsq)
                    ols_i_80_slopes.append(ols_i_slope)
                    ols_i_80_nu.append(red_nu)
                    ols_i_80_std.append(red_std)

def plot_slope_histograms():
    sns.set_style('white')
    fig, axs = plt.subplots(2,2)

    axs[0,0].hist(ols_68_slopes,color='cornflowerblue')
    axs[0,0].set_yticks([0,1,2])
    axs[0,0].set_xlabel('OLS Slope ('+r'$r^2$'+'>0.68)')
    axs[0,0].set_ylabel('Frequency')

    axs[0,1].hist(ols_80_slopes,color='cornflowerblue')
    axs[0,1].set_yticks([0,1,2])
    axs[0,1].set_xlabel('OLS Slope ('+r'$r^2$'+'>0.8)')

    axs[1,0].hist(ols_i_68_slopes,color='cornflowerblue')
    axs[1,0].set_xlabel('OLS Inliers Slope ('+r'$r^2$'+'>0.68)')
    axs[1,0].set_ylabel('Frequency')

    axs[1,1].hist(ols_i_80_slopes,color='cornflowerblue')
    axs[1,1].set_xlabel('OLS Inliers Slope ('+r'$r^2$'+'>0.8)')


    plt.subplots_adjust(wspace=0.25,hspace=0.4)
    plt.show()

def plot_comparisons():
    #OLS plot
    sns.set_style('white')

    fig, axs = plt.subplots(2,2)

    axs[0,0].scatter(ols_68_slopes,ols_68_nu,color='cornflowerblue')
    axs[0,0].set_xlabel('OLS Slope ('+r'$r^2$'+'>0.68)')
    axs[0,0].set_ylabel('Sokolovsky '+r'$\nu$')

    axs[0,1].scatter(ols_68_slopes,ols_68_std,color='cornflowerblue')
    axs[0,1].set_xlabel('OLS Slope ('+r'$r^2$'+'>0.68)')
    axs[0,1].set_ylabel(r'$\sigma$')

    axs[1,0].scatter(ols_80_slopes,ols_80_nu,color='cornflowerblue')
    axs[1,0].set_xlabel('OLS Slope ('+r'$r^2$'+'>0.80)')
    axs[1,0].set_ylabel('Sokolovsky '+r'$\nu$')

    axs[1,1].scatter(ols_80_slopes,ols_80_std,color='cornflowerblue')
    axs[1,1].set_xlabel('OLS Slope ('+r'$r^2$'+'>0.80)')
    axs[1,1].set_ylabel(r'$\sigma$')

    plt.subplots_adjust(wspace=0.4,hspace=0.4)
    plt.savefig('/home/thaddaeus/FMU/HRL/LAH/fifth_phase/for_email/summary_plots/ols_slope_comparisons.png')
    plt.clf()

    #OLS on inliers plot
    sns.set_style('white')

    fig, axs = plt.subplots(2,2)

    axs[0,0].scatter(ols_i_68_slopes,ols_i_68_nu,color='cornflowerblue')
    axs[0,0].set_xlabel('OLS Inlier Slope ('+r'$r^2$'+'>0.68)')
    axs[0,0].set_ylabel('Sokolovsky '+r'$\nu$')

    axs[0,1].scatter(ols_i_68_slopes,ols_i_68_std,color='cornflowerblue')
    axs[0,1].set_xlabel('OLS Inlier Slope ('+r'$r^2$'+'>0.68)')
    axs[0,1].set_ylabel(r'$\sigma$')

    axs[1,0].scatter(ols_i_80_slopes,ols_i_80_nu,color='cornflowerblue')
    axs[1,0].set_xlabel('OLS Inlier Slope ('+r'$r^2$'+'>0.80)')
    axs[1,0].set_ylabel('Sokolovsky '+r'$\nu$')

    axs[1,1].scatter(ols_i_80_slopes,ols_i_80_std,color='cornflowerblue')
    axs[1,1].set_xlabel('OLS Inlier Slope ('+r'$r^2$'+'>0.80)')
    axs[1,1].set_ylabel(r'$\sigma$')

    plt.subplots_adjust(wspace=0.4,hspace=0.4)
    plt.savefig('/home/thaddaeus/FMU/HRL/LAH/fifth_phase/for_email/summary_plots/ols_inliers_slope_comparisons.png')
    plt.clf()



plot_comparisons()















