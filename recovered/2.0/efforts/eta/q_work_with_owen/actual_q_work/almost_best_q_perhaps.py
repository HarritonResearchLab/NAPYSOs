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

    # Periodogram of the window function
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
    
    pwff = freq_W[np.where(power_W>high_power_W)] # pwff = Peak Window Function Frequencies
    
    wffitr = np.array([]) # wfftr = Window Function Frequency Indices To Remove

    ls_pers = 1/ls_freqs

    for f in pwff:
        per = 1/f
        #print((per-0.02),per,(per+0.02))
        #good_idx = np.invert(np.logical_and((ls_freqs<(f+dv)), (ls_freqs>(f-dv)))) # Phew that is clean
        good_idx = np.invert(np.logical_and(ls_pers<per+dv,ls_pers>per-dv))
        bap = np.where(good_idx==False)
        wffitr = np.append(wffitr,bap)#.astype(int)
        wffitr = np.unique(wffitr).astype(int)

        #print(wffitr)
    

    #print(wffitr)

    
    cleaned_powers = np.delete(ls_powers,wffitr)
    cleaned_freqs = np.delete(ls_freqs,wffitr)

    '''
    cleaned_powers = ls_powers[good_idx]
    cleaned_freqs = ls_freqs[good_idx]
    '''

    # Calculate FAPs
    faps = periodogram.false_alarm_level(false_alarm_levels)
    
    # Mask known aliase ranges (Lunar, etc.)
    cleaned_periods = 1/cleaned_freqs
    mask_lunar = np.invert(np.logical_and(cleaned_periods>26, cleaned_periods<30))
    cleaned_freqs = np.array(cleaned_freqs)[mask_lunar]
    cleaned_powers = np.array(cleaned_powers)[mask_lunar]

    # Find best results
    best_index = np.argmax(cleaned_powers)
    best_power = cleaned_powers[best_index]
    best_freq = cleaned_freqs[best_index]

    # Fold the light curve
    T = 1 / (float(best_freq))
    phased_dates = np.mod(mjd, T) / T  # Got this from the feets package documentation
    phased_dates_cycle_2 = phased_dates + 1

    def create_plot():
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
        ax1.set_xlim(0.5,250)
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
        create_plot()
        plt.show()
        plt.clf()
    
    elif out_type == 'None' or out_type == 'none':
        plt.clf()
    
    else:
        create_plot()
        plt.savefig(out_type, dpi=200, format=out_type[-3:]) # Flexible save type (svg, png, etc.)
        plt.clf()
    
    out_list = [best_power, 1/best_freq]
    for fap in faps: 
        out_list.append(fap)
    
    return out_list

def quas_per(mjd, mag, magerr, per, sig_factor):
    '''
    Function for calculating quasi-periodicity metric, originally adopted 
    from and heavily modified from Bredall et al. 2020's code. 
    '''
    # Import(s)
    import numpy as np
    import pandas as pd
    from astropy.convolution import Box1DKernel, convolve

    # Action
    JD = mjd + 2400000.5 

    # Calculate sig (!!!)
    sig = sig_factor*np.mean(magerr)

    # Create the residual curve
    phase = JD % per
    mag = mag[np.argsort(phase)]

    # We use three periods and extract the middle to prevent edge effects
    three_periods = np.concatenate((mag, mag, mag))
    boxcar = Box1DKernel(len(mag) // 4)
    smooth_mag = convolve(three_periods, boxcar)
    smooth_mag = smooth_mag[np.size(mag):2*np.size(mag)]

    resid_mag = mag - smooth_mag

    quas_per = ((np.nanstd(resid_mag) ** 2 - sig ** 2)/(np.nanstd(mag) ** 2 - sig ** 2))

    return quas_per, resid_mag

def codyM(x):
    #Import(s)
    from scipy import stats 
    import numpy as np
    
    #Action
    x = np.array(x)
    m_metric = (np.mean([stats.mstats.mquantiles(x,prob=0.9),stats.mstats.mquantiles(x,prob=0.1)])-np.median(x))/np.sqrt(((x-x.mean())**2).sum()/len(x))
    return m_metric

def test_new_q(key_file,data_path,plot_path,report_file):
    # Import(s)
    import pandas as pd
    import numpy as np
    from progress.bar import Bar

    # Action
    key_df = pd.read_csv(key_file)
    ids = list(key_df['ID'])

    out_ids = []
    qs = []
    ms = []
    
    bar = Bar('Processing...',max=(len(ids)-1))

    for id in ids:
        data_df = pd.read_csv(data_path.replace('+++',id))
        mjds = np.array(data_df['mjd'])
        mags = np.array(data_df['mag'])
        magerrs = np.array(data_df['magerr'])
        ls_res = lombScargle(mjd=mjds,mag=mags,magerr=magerrs,dv=0.25,min_per=0.5,max_per=250,false_alarm_levels=[0.1,0.05,0.01],out_type=plot_path)
        q = quas_per(mjd=mjds,mag=mags,magerr=magerrs,per=ls_res[1])[0]
        m = codyM(x=mags)
        qs.append(q)
        ms.append(m)
        out_ids.append(id)
        bar.next()
    
    bar.finish()
    zipped_list = list(zip(out_ids,qs,ms))
    out_df = pd.DataFrame(data=zipped_list,columns=['ID','Q','M'])   
    out_df.to_csv(report_file,index=False)

def q_diagnostic(id,mjd,mag,magerr,dv,min_per,max_per,false_alarm_levels,out_type):
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

    # Periodogram of the window function
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
    
    pwff = freq_W[np.where(power_W>high_power_W)] # pwff = Peak Window Function Frequencies
    
    wffitr = np.array([]) # wfftr = Window Function Frequency Indices To Remove

    ls_pers = 1/ls_freqs

    for f in pwff:
        per = 1/f
        #print((per-0.02),per,(per+0.02))
        #good_idx = np.invert(np.logical_and((ls_freqs<(f+dv)), (ls_freqs>(f-dv)))) # Phew that is clean
        good_idx = np.invert(np.logical_and(ls_pers<per+dv,ls_pers>per-dv))
        bap = np.where(good_idx==False)
        wffitr = np.append(wffitr,bap)#.astype(int)
        wffitr = np.unique(wffitr).astype(int)    
    
    cleaned_powers = np.delete(ls_powers,wffitr)
    cleaned_freqs = np.delete(ls_freqs,wffitr)

    '''
    cleaned_powers = ls_powers[good_idx]
    cleaned_freqs = ls_freqs[good_idx]
    '''

    # Calculate FAPs
    faps = periodogram.false_alarm_level(false_alarm_levels)
    
    # Mask known aliase ranges (Lunar, etc.)
    cleaned_periods = 1/cleaned_freqs
    mask_lunar = np.invert(np.logical_and(cleaned_periods>26, cleaned_periods<30))
    cleaned_freqs = np.array(cleaned_freqs)[mask_lunar]
    cleaned_powers = np.array(cleaned_powers)[mask_lunar]

    # Find best results
    best_index = np.argmax(cleaned_powers)
    best_power = cleaned_powers[best_index]
    best_freq = cleaned_freqs[best_index]
    best_per = 1/best_freq

    # Find second period and its power
    second_index = np.where(cleaned_powers==np.partition(cleaned_powers,-2)[2])
    second_power = float(cleaned_powers[second_index])
    second_freq = float(cleaned_freqs[second_index])
    second_per = 1/second_freq

    # Fold the light curve
    T = 1 / (float(best_freq))
    phased_dates = np.mod(mjd, T) / T  # Got this from the feets package documentation
    phased_dates_cycle_2 = phased_dates + 1

    # Calculate Q and get residuals plot
    qp_results = quas_per(mjd=mjd,mag=mag,magerr=magerr,per=best_per,sig_factor=1.25)
    q = qp_results[0]
    residuals = qp_results[1]

    # Calculate m
    m = codyM(x=mag)

    def create_plot():
        # Plot
        #mpl.style.use('classic')
        plt.rcParams['font.family'] = 'serif'
        fig, axs = plt.subplots(2,2)
        fig.suptitle('Q: '+str(round(q,3))+'; M: '+str(round(m,3)),fontsize='medium')

        # Periodogram plot
        periods = 1 / cleaned_freqs
        axs[1,1].plot(periods[np.argwhere(periods<26)], cleaned_powers[np.argwhere(periods<26)])
        axs[1,1].plot(periods[np.argwhere(periods>30)], cleaned_powers[np.argwhere(periods>30)],color='C0')
        
        #Plot FAP Levels
        colors = ['lightgrey','silver','darkgray','gray','dimgray']

        for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
            confidence_label = str(100*(1-false_alarm_levels[i]))[:-2]+'% FAP'
            axs[1,1].hlines(y=(faps[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
                label=confidence_label)

        axs[1,1].set_xscale('log')
        axs[1,1].set_xlabel('Period d')
        axs[1,1].set_xscale('log')
        axs[1,1].set_ylabel('Power')
        axs[1,1].set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
        axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())
        axs[1,1].set_ylim(0,1)
        axs[1,1].set_xlim(0.5,250)
        axs[1,1].set_title('Periodogram', fontsize=10)
        axs[1,1].legend(loc='upper right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
        axs[1,1].tick_params(axis='both',which='both',direction='in')
        #box = axs[1,1].get_position()
        #axs[1,1].set_position([box.x0,box.y0,box.width*0.5,box.height])
        #axs[1,1].legend(bbox_to_anchor=(1.15,0.5),loc='center',fontsize=4)

        # Folded light curve plot
        xlabel = 'Phase (P = '+str(round((1 / best_freq), 3))+' d)'
        axs[0,1].set_xlabel(xlabel)
        axs[0,1].set_ylabel('Mag')
        axs[0,1].scatter(phased_dates, mag, s=2)
        axs[0,1].scatter(phased_dates_cycle_2, mag, s=2, c='C0')
        axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
        axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
        axs[0,1].invert_yaxis()
        axs[0,1].set_title("Folded Light Curve", fontsize=10)
        axs[0,1].locator_params(axis='y', nbins=5)
        axs[0,1].tick_params(axis='both',which='both',direction='in')

        # Unfolded light curve plot
        axs[0,0].errorbar(mjd, mag, yerr=magerr, lw=0,elinewidth=0.5)
        axs[0,0].scatter(mjd,mag,s=2)
        axs[0,0].invert_yaxis()
        axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
        axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
        axs[0,0].set_ylabel('Mag')
        axs[0,0].set_xlabel('MJD')
        axs[0,0].set_title('Light Curve', fontsize=10)
        axs[0,0].locator_params(axis='y', nbins=5)
        axs[0,0].tick_params(axis='both',which='both',direction='in')

        # Residuals plot
        axs[1,0].scatter(mjd,residuals,s=2)
        axs[1,0].axhline(y=0,xmin=0,xmax=1,color='black',lw=0.75)
        axs[1,0].set_title('Residual Plot')
        axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
        axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())
        axs[1,0].set_xlabel('MJD')
        axs[1,0].set_ylabel('Residual Mag',fontsize=10)
        axs[1,0].tick_params(axis='both',which='both',direction='in')

        plt.subplots_adjust(wspace=0.4, hspace=0.50)

    if out_type == 'show' or out_type == 'Show':
        create_plot()
        plt.show()
        plt.clf()
    
    elif out_type == 'None' or out_type == 'none':
        plt.clf()
        plt.close()
    
    else:
        create_plot()
        q_str = str(q) 
        q_str_split = q_str.split('.')
        decimal_portion = q_str_split[1]
        decimal_portion = decimal_portion[0:3]
        q_str = q_str_split[0]+'.'+decimal_portion

        actual_path = out_type.replace('***',q_str)
        actual_path = actual_path.replace('+++',id)
        plt.savefig(actual_path, dpi=300, format=out_type[-3:]) # Flexible save type (svg, png, etc.)
        plt.clf()
        plt.close()
    
    out_list = [best_power, best_per,second_power,second_per,q,m]
    for fap in faps: 
        out_list.append(fap)
    
    return out_list

def run_q_diagnostic(key_file,data_path,plot_path,report_file):
    '''
    Create the q_diagnostic plots for all the objects. 
    '''
    # Import(s)
    import pandas as pd
    import numpy as np
    from progress.bar import Bar
    
    # Action
    key_df = pd.read_csv(key_file)
    ids = list(key_df['ID'])

    bar = Bar('Processing...',max=len(ids))
    with open(report_file,'a') as f:
        f.write('ID,POW1,PER1,POW2,PER2,Q,M,90%_FAP,95%_FAP,99%_FAP'+'\n')
    for id_index, id in enumerate(ids):
        if id_index>200:
            data_df = pd.read_csv(data_path.replace('+++',id))
            mjds = np.array(data_df['mjd'])
            mags = np.array(data_df['mag'])
            magerrs = np.array(data_df['magerr'])
            q_and_ls_res = q_diagnostic(id=id,mjd=mjds,mag=mags,magerr=magerrs,dv=0.25,min_per=0.5,max_per=250,false_alarm_levels=[0.1,0.05,0.01],out_type=plot_path)
            
            # Get results
            best_period_power = q_and_ls_res[0]
            best_period = q_and_ls_res[1]
            second_period_power = q_and_ls_res[2]
            second_period = q_and_ls_res[3]
            fitted_q = q_and_ls_res[4]
            fitted_m = q_and_ls_res[5]
            first_fap = q_and_ls_res[6]
            second_fap = q_and_ls_res[7]
            third_fap = q_and_ls_res[8]

            # Write results to line in report file
            result_list = [
                id,best_period_power,best_period,
                second_period_power,second_period,fitted_q,
                fitted_m,first_fap,second_fap,third_fap
                ]

            result_list = [str(val) for val in result_list]
            result_string = ','.join(result_list).replace('[','').replace(']','')
            
            with open(report_file,'a') as f:
                f.write(result_string+'\n')

        # Advance progress bar        
        bar.next()

    bar.finish()

dp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++_r.csv'
key_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/good_ids.csv'
plot_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/q_and_per_work/q_work/good_ids/plots/***:ID:+++.png'
res_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/q_and_per_work/q_work/good_ids/results.csv'
run_q_diagnostic(key_file=key_path,data_path=dp,plot_path=plot_temp,report_file=res_file)

