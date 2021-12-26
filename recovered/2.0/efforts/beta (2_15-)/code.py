"""
Note: JD = mjd + 2400000.5 
"""

def lombScargle(mjd,mag,magerr,dv,min_per,max_per,false_alarm_levels,out_type):
    #Import(s)
    from astropy.timeseries import LombScargle
    from pandas import Series, concat, date_range, to_datetime
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator

    #Action

    JD = mjd + 2400000.5

    # Generate a Dirac Comb, our window function
    time = to_datetime(JD, unit="D", origin="julian")
    time_not_obs = date_range(time.min(), time.max(), periods=1000)
    base = Series(np.zeros(len(time_not_obs)), index=time_not_obs)
    teeth = Series(np.ones(len(time)), index=time)
    dirac_comb = concat([base, teeth]).sort_index()
    
    minf = 1/max_per
    maxf = 1/min_per

    # First, the periodogram of the window function
    JD_W = dirac_comb.index.to_julian_date()
    mag_W = dirac_comb.values
    periodogram_W = LombScargle(JD_W, mag_W)
    freq_W, power_W = periodogram_W.autopower(method='fastchi2',minimum_frequency=minf, maximum_frequency=maxf)


    # Periodogram of original light curve
    periodogram = LombScargle(JD, mag, magerr)
    ls_freqs, ls_powers = periodogram.autopower(method='fastchi2',minimum_frequency=minf, maximum_frequency=maxf)

    # Mask out peak window-function frequencies from the data with a notch
    # width of dv (default should be 0.03) Hz on either side.
    high_power_W = power_W.mean() + 2 * power_W.std()
    

    pwf = freq_W[np.argwhere(power_W>high_power_W)] # pwf = Peak Window Frequencies
    indices_to_remove = []

    for f in pwf: 
        mask = np.invert(np.logical_and((ls_freqs+dv)>f, (ls_freqs-dv)<f))
        idx = np.argwhere(mask==False)
        for i in idx:
            indices_to_remove.append(int(i))
    
    indices_to_remove = np.array(sorted(indices_to_remove,reverse=True))


    cleaned_freqs = np.delete(ls_freqs,indices_to_remove)
    cleaned_powers = np.delete(ls_powers,indices_to_remove)
    
    # Calculate FAPs
    faps = periodogram.false_alarm_level(false_alarm_levels)
    
    cleaned_periods = 1/cleaned_freqs
    mask_lunar = np.invert(np.logical_and(cleaned_periods>26, cleaned_periods<30))
    cleaned_freqs = cleaned_freqs[mask_lunar]
    cleaned_powers = cleaned_powers[mask_lunar]
    cleaned_periods = cleaned_periods[mask_lunar]

    '''
    day_mask_one = np.invert(np.logical_and(cleaned_periods>0.5, cleaned_periods<0.505)) # First "day" mask, exclude 0.5<P<0.505
    cleaned_freqs = cleaned_freqs[day_mask_one]
    cleaned_powers = cleaned_powers[day_mask_one]
    cleaned_periods = cleaned_periods[day_mask_one]

    day_mask_two = np.invert(np.logical_and(cleaned_periods>0.97, cleaned_periods<1.04)) # Second "day" mask, exclude 0.97<P<1.04
    cleaned_freqs = cleaned_freqs[day_mask_two]
    cleaned_powers = cleaned_powers[day_mask_two]
    cleaned_periods = cleaned_periods[day_mask_two]
    '''
    min_mask = np.logical_and(cleaned_periods>1.5,cleaned_periods<=250)
    cleaned_freqs = cleaned_freqs[min_mask]
    cleaned_powers = cleaned_powers[min_mask]
    cleaned_periods = cleaned_periods[min_mask]

    best_index = np.argmax(cleaned_powers)
    best_power = cleaned_powers[best_index]
    best_freq = cleaned_freqs[best_index]

    # Fold the light curve
    T = 1 / (float(best_freq))
    phased_dates = np.mod(mjd, T) / T  # Got this from the feets package documentation
    phased_dates_cycle_2 = phased_dates + 1

    # Plot
    plt.rcParams['font.family'] = 'serif'
    ax1 = plt.subplot(222)

    #first_plot_mask = np.logical_and(periods>0.505, periods<0.97)
    #second_plot_mask = np.logical_and(periods>1.04,periods<26)
    
    third_plot_mask = np.argwhere(cleaned_periods>30)
    fourth_plot_mask = np.argwhere(cleaned_periods<26)

    #ax1.plot(periods[first_plot_mask], cleaned_powers[first_plot_mask])
    #ax1.plot(periods[second_plot_mask], cleaned_powers[second_plot_mask],color='C0')
    ax1.plot(cleaned_periods[third_plot_mask], cleaned_powers[third_plot_mask],color='C0')
    ax1.plot(cleaned_periods[fourth_plot_mask], cleaned_powers[fourth_plot_mask],color='C0')
    
    #Plot FAP Levels
    colors = ['lightgrey','silver','darkgray','gray','dimgray']

    for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
        confidence_label = str(100*(1-false_alarm_levels[i]))+'% FAP'
        ax1.hlines(y=(faps[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
               label=confidence_label)

    # Periodogram plot 
    ax1.set_xscale('log')
    ax1.set_xlabel('Period d')
    ax1.set_xscale('log')
    ax1.set_ylabel('Power')
    ax1.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.set_ylim(0,1)
    ax1.set_xlim(1.5,250)
    ax1.set_title('Periodogram', fontsize=10)
    #box = ax1.get_position()
    #ax1.set_position([box.x0,box.y0,box.width*0.5,box.height])
    #ax1.legend(bbox_to_anchor=(1.15,0.5),loc='center',fontsize=4)
    ax1.legend(loc='upper right',fontsize=4)

    # Folded light curve plot
    ax2 = plt.subplot(221)
    xlabel = 'Phase (P = '+str(round((1/best_freq), 3))+' d)'
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel('Mag')
    ax2.scatter(phased_dates, mag, s=2)
    ax2.scatter(phased_dates_cycle_2, mag, s=2, c='C0')
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.invert_yaxis()
    ax2.set_title("Folded Light Curve", fontsize=10)
    ax2.locator_params(axis='y', nbins=5)

    # Unfolded light curve plot
    ax3 = plt.subplot(212)
    ax3.errorbar(mjd, mag, yerr=magerr, lw=0,elinewidth=0.5)
    ax3.scatter(mjd,mag,s=2)
    ax3.invert_yaxis()
    ax3.xaxis.set_minor_locator(AutoMinorLocator())
    ax3.yaxis.set_minor_locator(AutoMinorLocator())
    ax3.set_ylabel('Mag')
    ax3.set_xlabel('MJD')
    ax3.set_title('Light Curve', fontsize=10)
    ax3.locator_params(axis='y', nbins=5)

    plt.subplots_adjust(wspace=0.4, hspace=0.45)

    if out_type == 'show' or out_type == 'Show':
        plt.show()
        plt.clf()
    else:
        plt.savefig(out_type, dpi=200, format=out_type[-3:]) # Flexible save type (svg, png, etc.)
        plt.clf()
    
    out_list = [1/best_freq, best_power]
    for fap in faps: 
        out_list.append(fap)
    
    return out_list
    
def test(data_file):
    #Import(s)
    import pandas as pd
    import numpy as np

    #Action
    df = pd.read_csv(data_file)
    mjds = np.array(df['mjd'])
    mags = np.array(df['mag'])
    magerrs = np.array(df['magerr'])

    bp2 = lombScargle(mjd=mjds,mag=mags,magerr=magerrs,dv=0.03,min_per=0.5,max_per=250,false_alarm_levels=[0.1,0.05,0.01],out_type='Show')

    print(bp2[0],bp2[1],bp2[2])

def redo_lomb_on_all_again(key,data_path,plot_path,out_path):
    #Import(s)
    import pandas as pd
    import numpy as np
    from progress.bar import ShadyBar

    #Action
    key_df = pd.read_csv(key)
    ids = list(key_df['ID'])
    bar = ShadyBar('Processing...',max=len(ids))

    out_ids = []
    best_pers = []
    best_pows = []
    fap1 = []
    fap2 = []
    fap3 = []
    fap4 = []

    for id in ids:
        data_df = pd.read_csv(data_path.replace('+++',id))
        mjds = np.array(data_df['mjd'])
        mags = np.array(data_df['mag'])
        magerrs = np.array(data_df['magerr'])

        save_path = plot_path.replace('+++',id)

        ls = lombScargle(mjd=mjds,mag=mags,magerr=magerrs,dv=0.03,min_per=0.5,max_per=250,false_alarm_levels=[0.1,0.05,0.01,0.001],out_type=save_path)
        if ls[0] != np.nan:
            out_ids.append(id)
            best_pers.append(ls[0])
            best_pows.append(ls[1])
            fap1.append(ls[2])
            fap2.append(ls[3])
            fap3.append(ls[4])
            fap4.append(ls[5])
        
        bar.next()
    
    zipped = list(zip(out_ids,best_pers,best_pows,fap1,fap2,fap3,fap4))
    out_df = pd.DataFrame(data=zipped,columns=['ID','PER','POW','90FAP','95FAP','99FAP','999FAP'])
    out_df.to_csv(out_path)

    bar.finish()

key_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/key.csv'
r_path_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++_r.csv'
plot_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/delta/stewardship_stuff/redone_ls_again/plots/+++.svg' 
report_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/delta/stewardship_stuff/redone_ls_again/report_file.csv'

redo_lomb_on_all_again(key=key_file,data_path=r_path_temp,plot_path=plot_path,out_path=report_file)
test('/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/FHK_118_r.csv')