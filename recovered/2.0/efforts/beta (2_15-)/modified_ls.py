def lombScargle(x, y, yerr, fap_levels, min_period, out_type):
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



def redo_lomb_on_all():
    #Import(s)
    import pandas as pd
    from progress.bar import Bar

    #Action
    key = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/key.csv'

    df = pd.read_csv(key)
    ids = list(df['ID'])
    periods = []
    powers = []
    fap_0 = []
    fap_1 = []
    fap_2 = []
    fap_3 = []

    bar = Bar('Processing...',max=len(ids))
    for id in ids:
        r_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++_r.csv'.replace('+++',id)
        data_df = pd.read_csv(r_file)
        r_mags = list(data_df['mag'])
        r_mjds = list(data_df['mjd'])
        r_magerrs = list(data_df['magerr'])

        plot_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/gamma/updating_lomb_scarg/plots/+++.svg'.replace('+++',id)

        ls = lombScargle(x=r_mjds,y=r_mags,yerr=r_magerrs,fap_levels=[0.1, 0.05, 0.01, 0.001],min_period=1.81,out_type=plot_path)
        periods.append(ls[0])
        powers.append(ls[1])
        fap_0.append(ls[2])
        fap_1.append(ls[3])
        fap_2.append(ls[4])
        fap_3.append(ls[5])
        bar.next()
    
    zipped_list = list(zip(ids,periods,powers,fap_0,fap_1,fap_2,fap_3))
    save_df = pd.DataFrame(data=zipped_list,columns=['ID','Period','Power','90% FAP','95% FAP','99% FAP','99.9% FAP'])
    save_df.to_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/gamma/updating_lomb_scarg/updated_ls_results.csv',index=False)

redo_lomb_on_all()

