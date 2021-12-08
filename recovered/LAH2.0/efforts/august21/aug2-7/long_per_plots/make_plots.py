# Import(s)
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Action
 
def make_plot(id, plot_dir, mjd, mag, magerr): 
    # Calculate periodogram
    fap_levels = [0.1, 0.05, 0.01]
    frequencies, powers = LombScargle(mjd, mag, magerr).autopower(minimum_frequency=(1 / 250),
                                                             maximum_frequency=(1 / 0.55))
    # Some statistics
    false_alarm_levels = list(LombScargle(mjd, mag, magerr).false_alarm_level(fap_levels))
    
    # Phase the data    
    best_frequency = frequencies[np.argmax(powers)]
    T = 1 / (float(best_frequency))
    phased_dates = np.mod(mjd, T) / T  # grabbed this from feets package documentation
    phased_dates_cycle_2 = phased_dates + 1
    
    periods = 1 / frequencies
    
    # get long period
    periods_mask = np.where(periods>25)
    periods_test = periods[periods_mask]
    powers_test = powers[periods_mask]
    long_per = periods_test[np.argsort(powers_test)[-1]]
    
    # Actually make plot 
    
    plt.rcParams['font.family']='serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    plt.rcParams['font.size']=9
    
    # bottom long
    ax1 = plt.subplot(212)
    ax1.scatter(mjd, mag, color='#408ee0', edgecolors='black', s=6, lw=0.25)
    ax1.set_ylabel('Mag')
    ax1.set_xlabel('MJD')
    ax1.invert_yaxis()
    ax1.set_title('Light Curve', fontsize=9)
    
    # Top left
    ax2 = plt.subplot(231)
    
    xlabel = 'Phase (P = '+str(round((1 / best_frequency), 3))+' d)'
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel('Mag')
    ax2.scatter(phased_dates, mag, color='#408ee0', edgecolors='black', s=6, lw=0.25)
    ax2.scatter(phased_dates_cycle_2, mag, color='#408ee0', edgecolors='black', s=6, lw=0.25)
    ax2.invert_yaxis()
    ax2.set_title('Folded Light Curve', fontsize=9)
    
    # Top middle
    ax3 = plt.subplot(232)
    ax3.plot(periods, powers, color='#408ee0')
    #Plot FAP Levels
    colors = ['lightgrey','silver','darkgray','gray','dimgray']

    for i in range(len(fap_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
        confidence_label = str(100*(1-fap_levels[i]))+'% FAP'
        ax3.hlines(y=(false_alarm_levels[i]), xmin=0.55, xmax=250, color = colors[i],linestyles='--',
               label=confidence_label)

    ax3.set_xscale('log')
    ax3.set_xlabel('Period d')
    ax3.set_xscale('log')
    ax3.set_ylabel('Power')
    #ax1.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
    #ax3.margins(.1, .1)
    ax3.yaxis.set_minor_locator(AutoMinorLocator())
    #ax3.legend(fontsize=6, loc='upper center', edgecolor='black')
    ax3.set_title('Periodogram', fontsize=9)
    
    # Top right
    ax4 = plt.subplot(233)
    phased_dates = np.mod(mjd, long_per) / long_per  # grabbed this from feets package documentation
    phased_dates_cycle_2 = phased_dates + 1
    xlabel = 'Phase (P = '+str(round(long_per, 3))+' d)'
    ax4.set_xlabel(xlabel)
    ax4.set_ylabel('Mag')
    ax4.scatter(phased_dates, mag, color='#408ee0', edgecolors='black', s=6, lw=0.25)
    ax4.scatter(phased_dates_cycle_2, mag, color='#408ee0', edgecolors='black', s=6, lw=0.25)
    ax4.invert_yaxis()
    ax4.set_title('Folded Light Curve', fontsize=9)
    
    # Shared settings
    for ax in [ax1, ax2, ax4]:
        ax.tick_params(axis='both',which='both',direction='in')
        ax.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        ax.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    fig = plt.gcf()
    fig.set_size_inches(8.4, 5.2)
    
    plt.subplots_adjust(wspace=0.5, hspace=0.45)
    
    if plot_dir == 'show':
        plt.show()
    else: 
        plt.savefig(plot_dir+'/'+id+'.png', dpi=300)
    plt.clf()
    plt.close()

first_df = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/light_curves/FHK_286_r.csv')
first_mjds = np.array(first_df['mjd'])
first_mags = np.array(first_df['mag'])
first_magerrs = np.array(first_df['magerr'])

plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/aug2-7/long_per_plots/plots'
make_plot('FHK_286', plot_dir, first_mjds, first_mags, first_magerrs)

second_df = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/FHK_176/FHK_176_r.csv')
second_mjds = np.array(second_df['mjd'])
second_mags = np.array(second_df['mag'])
second_magerrs = np.array(second_df['magerr'])
make_plot('FHK_176', plot_dir, second_mjds, second_mags, second_magerrs)