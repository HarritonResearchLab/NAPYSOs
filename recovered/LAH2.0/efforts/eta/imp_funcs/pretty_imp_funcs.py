def fig_plot_2(data_file,q_bounds):
    '''
    Create figure two plot (from Bredall et al. 2020) from data file which has all the IDs, Qs, and Ms 
    '''
    
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.ticker import MultipleLocator

    # Action

    df = pd.read_csv(data_file)
    ids = list(df['ID'])
    q_raw = np.array(df['Q'])
    for i in q_raw:
        print(i)
    q_mask = np.logical_and(q_raw>0.0,q_raw<1.0)
    
    m = np.array(df['M'])[q_mask]
    q_raw = q_raw[q_mask]
    
    # Transform Q
    q = (q_raw-min(q_raw))*(1/(max(q_raw)-min(q_raw)))
    lower_q_bound = (q_bounds[0]-min(q_raw))*(1/(max(q_raw)-min(q_raw)))
    upper_q_bound = (q_bounds[1]-min(q_raw))*(1/(max(q_raw)-min(q_raw)))
    #q = q_raw

    # Get morph classes 

    dipper_mask = np.where(m>0.25)
    dipper_ms = m[dipper_mask]
    dipper_qs = q[dipper_mask]

    burster_mask = np.where(m<-0.25)
    burster_ms = m[burster_mask]
    burster_qs = q[burster_mask]

    sym_mask = np.logical_and(m<=0.25,m>=-0.25)
    sym_ms = m[sym_mask]
    sym_qs = q[sym_mask]

    
    # Plot
    plt.rcParams['font.family'] = 'serif'
    plt.gca().tick_params(axis='both',which='both',direction='in')
    plt.scatter(dipper_qs,dipper_ms,marker='o',color='#408ee0',s=4) 
    plt.scatter(sym_qs,sym_ms,marker='s',color='#89bff8',s=4)
    plt.scatter(burster_qs,burster_ms,marker='d',color='#0051a2',s=4)

    plt.figtext(0.18,0.9,'Periodic',ha='center',fontsize='small')
    plt.figtext(0.5,0.9,'Quasi-Periodic',ha='center',fontsize='small')
    plt.figtext(0.845,0.9,'Aperiodic',ha='center',fontsize='small')

    plt.figtext(0.91,0.685,'Bursting',rotation=270,fontsize='small')
    plt.figtext(0.91,0.43,'Symmetric',rotation=270,fontsize='small')
    plt.figtext(0.91,0.22,'Dipping',rotation=270,fontsize='small')

    plt.axhline(y=0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axhline(y=-0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=lower_q_bound,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=upper_q_bound,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.xlim(0,1)
    plt.ylim(-1,1)
    plt.gca().invert_yaxis()
    plt.xlabel('Quasi-Periodicity (Q)')
    plt.ylabel('Flux Asymmetry (M)')
    plt.axes().yaxis.set_minor_locator(MultipleLocator(0.05))
    plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
    plt.tick_params(axis='both',which='both',direction='in')
    plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)

    plt.show()

def lombScargle(mjd,mag,magerr,dv,min_per,max_per,false_alarm_levels,out_type):
    ''' 
    Best working version of our lomb-scargle function
    '''
    
    # Import(s)
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

    for f in pwf: 
        good_idx = np.invert(np.logical_and((ls_freqs+dv)>f, (ls_freqs-dv)<f)) # Phew that is clean
    
    cleaned_powers = ls_powers[good_idx]
    cleaned_freqs = ls_freqs[good_idx]
    
    # Calculate FAPs
    faps = periodogram.false_alarm_level(false_alarm_levels)
    
    cleaned_periods = 1/cleaned_freqs
    mask_lunar = np.invert(np.logical_and(cleaned_periods>26, cleaned_periods<30))
    cleaned_freqs = np.array(cleaned_freqs)[mask_lunar]
    cleaned_powers = np.array(cleaned_powers)[mask_lunar]
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
    periods = 1 / cleaned_freqs
    ax1.plot(periods[np.argwhere(periods<26)], cleaned_powers[np.argwhere(periods<26)])
    ax1.plot(periods[np.argwhere(periods>30)], cleaned_powers[np.argwhere(periods>30)],color='C0')
    
    #Plot FAP Levels
    colors = ['lightgrey','silver','darkgray','gray','dimgray']

    for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
        confidence_label = str(100*(1-false_alarm_levels[i]))+'% FAP'
        ax1.hlines(y=(false_alarm_levels[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
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
    box = ax1.get_position()
    ax1.set_position([box.x0,box.y0,box.width*0.5,box.height])
    ax1.legend(bbox_to_anchor=(1.15,0.5),loc='center',fontsize=4)

    # Folded light curve plot
    ax2 = plt.subplot(221)
    xlabel = 'Phase (P = '+str(round((1 / best_freq), 3))+' d)'
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
    
    out_list = [best_power, 1/best_freq]
    for fap in faps: 
        out_list.append(fap)
    
    return out_list

fig_plot_2(data_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results.csv',q_bounds=[0.45,0.8])