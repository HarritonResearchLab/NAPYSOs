def bredall_best_per(JD, mag, err):
    #Import(s)
    from astropy.timeseries import LombScargle
    from pandas import Series, concat, date_range, to_datetime
    import numpy as np

    #Action

    # Generate a Dirac Comb, our window function
    time = to_datetime(JD, unit="D", origin="julian")
    time_not_obs = date_range(time.min(), time.max(), periods=1000)
    base = Series(np.zeros(len(time_not_obs)), index=time_not_obs)
    teeth = Series(np.ones(len(time)), index=time)
    dirac_comb = concat([base, teeth]).sort_index()

    # We only want to look for signals that we have 3 periods observed. ###I SHOULD ADD THIS TO OUR LS FUNC
    maxp = (JD.max() - JD.min()) / 3
    minf = 1/maxp

    # First, we find the periodigram of the window function
    JD_W = dirac_comb.index.to_julian_date()
    mag_W = dirac_comb.values
    periodogram_W = LombScargle(JD_W, mag_W)
    freq_W, power_W = periodogram_W.autopower(minimum_frequency=minf, maximum_frequency=1/0.1)

    # Now, for the original lightcurve
    periodogram = LombScargle(JD, mag)
    freq, power = periodogram.autopower(minimum_frequency=minf, maximum_frequency=1/0.1)

    # Mask out peak window-function frequencies from the data with a notch
    # width of 0.03 Hz on either side.
    high_power_W = power_W.mean() + 2 * power_W.std()
    for idx in np.argwhere(power_W > high_power_W):
        v = freq[idx]
        dv = 0.03
        vmin, vmax = v - dv, v + dv
        power[np.logical_and(freq > vmin, freq < vmax)] = 0

    # We find the best frequency that has a false alarm probability < 0.1.
    faps = periodogram.false_alarm_probability(power)
    try: 
        idx = np.argmax(power)
        best_fap = faps[idx]
        while best_fap > 0.1:
            power = np.delete(power, idx)
            idx = np.argmax(power)
            best_fap = faps[idx]
        best_freq = freq[idx]

        # Periods are more human-friendly :)
        return 1 / best_freq
    except ValueError:
        return 'NaN'
    
def ourLombScargle(x, y, yerr, fap_levels, min_period, out_type):
    # Import(s)
    from astropy.timeseries import LombScargle
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib

    # Action
    dates = np.array(x)
    mags = np.array(y)

    # Calculate periodogram
    max_period = (max(x)-min(x))/3
    frequencies, powers = LombScargle(dates, mags, yerr).autopower(method='fastchi2', minimum_frequency=(1 / max_period),
                                                             maximum_frequency=(1 / min_period))
    # Some statistics
    false_alarm_levels = list(LombScargle(dates, mags, yerr).false_alarm_level(fap_levels))
    # Phase the data
    best_frequency = frequencies[np.argmax(powers)]
    best_power = powers[np.argmax(powers)]
    T = 1 / (float(best_frequency))
    phased_dates = np.mod(dates, T) / T  # grabbed this from feets package documentation
    phased_dates_cycle_2 = phased_dates + 1

    # Generate plots
    plt.rcParams['font.family'] = 'serif'
    ax1 = plt.subplot(222)
    periods = 1 / frequencies
    ax1.plot(periods, powers)
    
    #Plot FAP Levels
    colors = ['lightgrey','silver','darkgray','gray','dimgray']

    for i in range(len(fap_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
        confidence_label = str(100*(1-fap_levels[i]))+'% FAP'
        ax1.hlines(y=(false_alarm_levels[i]), xmin=min_period, xmax=max_period, color = colors[i],linestyles='--',
               label=confidence_label)

    ax1.set_xscale('log')
    ax1.set_xlabel('Period d')
    ax1.set_xscale('log')
    ax1.set_ylabel('Power')
    #ax1.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax1.set_title('Periodogram', fontsize=10)
    box = ax1.get_position()
    ax1.set_position([box.x0,box.y0,box.width*0.5,box.height])
    ax1.legend(bbox_to_anchor=(1.15,0.5),loc='center',fontsize=4)

    ax2 = plt.subplot(221)
    xlabel = 'Phase (P = '+str(round((1 / best_frequency), 3))+' d)'
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel('Mag')
    ax2.scatter(phased_dates, mags, s=2)
    ax2.scatter(phased_dates_cycle_2, mags, s=2, c='C0')
    ax2.invert_yaxis()
    ax2.set_title("Folded Light Curve", fontsize=10)
    ax2.locator_params(axis='y', nbins=5)

    ax3 = plt.subplot(212)
    ax3.errorbar(dates, mags, yerr=yerr, lw=0,elinewidth=0.5)
    ax3.scatter(dates,mags,s=2)
    ax3.invert_yaxis()
    ax3.set_ylabel('Mag')
    ax3.set_xlabel('MJD')
    ax3.set_title('Light Curve', fontsize=10)
    ax3.locator_params(axis='y', nbins=5)

    plt.subplots_adjust(wspace=0.4, hspace=0.45)

    if out_type == 'show':
        plt.show()
        plt.clf()
    else:
        plt.savefig(out_type, dpi=200, format='svg')
        plt.clf()
    
    out_list = [(1/best_frequency),best_power]
    for fap in false_alarm_levels:
        out_list.append(float(fap))

    return out_list