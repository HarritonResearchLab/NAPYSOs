def june8thQ(key,r_data_temp, g_data_temp, min_per,max_per,false_alarm_levels, out_type, report_file):
    ''' 
    ***
    Best working version of our lomb-scargle, flux-asymmetry, and quasi-periodicity code
    Certain sections were adopted and subsequently modified from code used in Bredall
    et al. (2020). If the dominant, aka peak r band period is between 0.5-0.51 or 0.99-1.01, 
    then this function removes periods in those ranges and re-evaluates what the dominant 
    period is. If no significant r band period is found, this function returns the dominant
    g band period (if it is significant). If no significant period is found, or if the 
    period is more than 249 days (e.g. 250.0), then the period is returned as np.nan.  
    ***
    key: .csv file with a column of object ids labeled "ID"  
    r_data_temp: path to g data file, with object id replaced with the string '+++' 
    g_data_temp: path to g data file, "" 
    min_per: minimum period. Set to 0.5
    max_per: maximum period. Set to 250
    false_alarm_levels: FAPs, see documentation for astropy Lomb-Scargle.
                        Set to [0.1, 0.05, 0.01] for 90%, 95%, and 99% 
                        confidence evaluations. 
    out_type: 'Show' shows the plot, 'none' or 'None' skips creating the plot, 
            and anything else is used as the directory into which the plot
            should be saved. Note that this should be the directory, not the usual 
            file path with the id replaced with '+++'. 
    report_file: the full path for the file into which the rest (majority)
                 of the results should be deposited into
    '''
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    from scipy.signal import find_peaks
    from progress.bar import Bar 

    # Some function definitions
    def codyM(x):
        #Import(s)
        from scipy import stats 
        import numpy as np
        
        #Action
        x = np.array(x)
        m_metric = (np.mean([stats.mstats.mquantiles(x,prob=0.9),stats.mstats.mquantiles(x,prob=0.1)])-np.median(x))/np.sqrt(((x-x.mean())**2).sum()/len(x))
        
        return m_metric

    def quas_per(mjd, mag, magerr, per, sig_factor):
        '''
        Function for calculating quasi-periodicity metric 
        '''
        # Import(s)
        import numpy as np
        from astropy.convolution import Box1DKernel, convolve

        # Action

        # Calculate sig
        sig = sig_factor*np.mean(magerr)

        # Create the residual curve
        phase = mjd % per
        mag = mag[np.argsort(phase)]

        # We use three periods and extract the middle to prevent edge effects
        three_periods = np.concatenate((mag, mag, mag))
        boxcar = Box1DKernel(len(mag) // 4)
        smooth_mag = convolve(three_periods, boxcar)
        smooth_mag = smooth_mag[np.size(mag):2*np.size(mag)]

        resid_mag = mag - smooth_mag

        quas_per = ((np.nanstd(resid_mag) ** 2 - sig ** 2)/(np.nanstd(mag) ** 2 - sig ** 2))

        return float(quas_per), resid_mag

    def sokolovskyNu(x,xerr):
        '''
        Best working version for calculating nu variabilty metric 
        from Sokolovsky et al. 2017. 
        '''
        
        #Import(s)
        import numpy as np
        
        #Action
        # [ (m-e)max - (m+e)min ] / [ (m-e)max + (m+e)min ]

        x = np.array(x)
        xerr = np.array(xerr)
        
        x_minus_xerr = x-xerr
        x_plus_xerr = x+xerr
        
        min_val = np.amin(x_plus_xerr)
        max_val = np.amax(x_minus_xerr)
        
        nu = (max_val-min_val)/(max_val+min_val)
        return nu
          
    def lomb_scargle_routine(mjd, mag, magerr, minf, maxf): 
        # Import(s)
        from astropy.timeseries import LombScargle
        
        # Action
        
        # Periodogram of light curve
        periodogram = LombScargle(mjd,mag,magerr)
        all_freqs, all_powers = periodogram.autopower(minimum_frequency=minf,maximum_frequency=maxf)
        faps = periodogram.false_alarm_level(false_alarm_levels)
        
        all_periods = 1/all_freqs
        
        # Find top five periods
        peak_power_indices, _ = find_peaks(all_powers, distance=5)
        all_peak_powers = np.sort(all_powers[peak_power_indices])[::-1]
        
        # Periods will be sorted by their powers, 
        # most powerful period coming first
        
        peak_powers = np.array([])
        peak_periods = np.array([])
                
        for i in range(0,5): 
            peak_index = np.where(all_powers==all_peak_powers[i])
            peak_powers = np.append(peak_powers, all_powers[peak_index])
            peak_periods = np.append(peak_periods, all_periods[peak_index][0])
             
        # Check for flag criteria 
        def flag(period): 
            if 0.99 < period < 1.01: 
                return True
            elif 0.5 < period < 0.51: 
                return True
            else: 
                return False
        
        return peak_powers, peak_periods, all_powers, all_periods, faps, flag(peak_periods[0]) 
     
    def return_cleaned(periods, powers): 
        # Cleaned as in the points in the bad ranges
        # e.g. 0.5-0.51, 0.99-1.01 have been deleted
        # Only used for objects that have peak periods in 
        # these intervals
        
        for range_to_ignore in ['0.5:0.51','0.99:1.01']:
            range_to_ignore = range_to_ignore.split(':')
            lower_bound = float(range_to_ignore[0])
            upper_bound = float(range_to_ignore[1])
            indices_to_delete = np.logical_and(periods<upper_bound,periods>lower_bound)
            periods = np.delete(arr=periods,obj=indices_to_delete)
            powers = np.delete(arr=powers,obj=indices_to_delete)
        
        return periods, powers
     
    def make_plot(mjds, mags, magerrs, peak_per, all_pers, all_powers, faps, Q, M):
        # Import(s)
        from matplotlib import ticker
        from matplotlib.ticker import AutoMinorLocator
        
        # Definitions
        def fold_lc(mjds, per):
            phased_dates = (np.mod(mjds, per))/per 
            phased_dates_cycle_2 = phased_dates + 1
            return phased_dates, phased_dates_cycle_2
        
        # Action 
        plt.rcParams['font.family'] = 'serif'
        plt.suptitle('Q: '+str(round(Q,3))+'; M: '+str(round(M, 3)),fontsize=10)
        
        ax1 = plt.subplot(221) # folded LC
        folded_dates = fold_lc(mjds, peak_per)
        ax1.scatter(folded_dates[0], mags, color='#408ee0', s=3)
        ax1.scatter(folded_dates[1], mags, color='#408ee0', s=3)
        ax1.set_xlabel('Phase (P='+str(round(peak_per,2))+' d)')
        ax1.set_ylabel('mag')
        ax1.set_title('Folded Light Curve', fontsize='medium')
        ax1.tick_params(axis='both',which='both',direction='in')
        ax1.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        ax1.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.yaxis.set_minor_locator(AutoMinorLocator())
        ax1.invert_yaxis()
        
        ax2 = plt.subplot(222) # Periodogram 
        ax2.plot(all_pers,all_powers,color='#408ee0')
        
        #Plot FAP Levels
        colors = ['lightgrey','silver','darkgray','gray','dimgray']

        for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
            confidence_label = str(100*(1-false_alarm_levels[i]))[:-2]+'% FAP'
            ax2.hlines(y=(faps[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
                label=confidence_label)

        ax2.set_xscale('log')
        ax2.set_xlabel('Period d')
        ax2.set_xscale('log')
        ax2.set_ylabel('Power')
        ax2.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
        ax2.set_ylim(0,1)
        ax2.set_xlim(min_per,max_per)
        ax2.set_title('Periodogram', fontsize='medium')
        ax2.tick_params(axis='both',which='both',direction='in')
        ax2.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        ax2.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        locmaj = ticker.LogLocator(base=10.0, subs=(0.158, 0.2511, 0.398,0.63 ))
        ax2.xaxis.set_minor_locator(locmaj)
        ax2.yaxis.set_minor_locator(AutoMinorLocator())
        ax2.legend(loc='upper right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
        
        ax3 = plt.subplot(212) # Unfolded LC
        ax3.errorbar(mjds, mags, yerr=magerrs, lw=0,elinewidth=0.5, color = "#408ee0")
        ax3.scatter(mjds,mags,s=2, color = "#408ee0")
        ax3.invert_yaxis()
        ax3.set_ylabel('Mag')
        ax3.set_xlabel('MJD')
        ax3.set_title('Light Curve', fontsize='medium')
        ax3.locator_params(axis='y', nbins=5)
        ax3.tick_params(axis='both',which='both',direction='in')
        ax3.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        ax3.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        ax3.xaxis.set_minor_locator(AutoMinorLocator())
        ax3.yaxis.set_minor_locator(AutoMinorLocator())
        
        plt.subplots_adjust(hspace=0.3, wspace=0.4)
    
    # Action

    key_df = pd.read_csv(key)
    ids = list(key_df['ID'])
    ras = np.array(key_df['RA'])
    decs = np.array(key_df['DEC'])
    
    bar = Bar('Processing...', max=len(ids))

    if report_file != 'None' or report_file != 'none':
        with open(report_file,'a') as f:
            column_names = 'ID,RA,DEC,PER,Q,M,NU,rMean'
            f.write(column_names+'\n')
           
    for id, ra, dec in zip(ids, ras, decs):
            # Run analysis for red data
            r_data_file = r_data_temp.replace('+++',id)
            r_data_df = pd.read_csv(r_data_file)
            r_mjd = np.array(r_data_df['mjd'])
            r_mag = np.array(r_data_df['mag'])
            r_magerr = np.array(r_data_df['magerr'])

            minf = 1/max_per
            maxf = 1/min_per 

            r_ls = lomb_scargle_routine(r_mjd, r_mag, r_magerr, minf, maxf)
            
            r_period_flag = r_ls[5]
            
            if r_period_flag == True: 
                r_faps = r_ls[4]
                
                r_cleaned_results = return_cleaned(r_ls[3], r_ls[2])
            
                r_all_powers = r_cleaned_results[1]
                r_all_periods = r_cleaned_results[0]
                
                r_sorted_powers = np.sort(r_all_powers)[::-1]
                r_peak_indices = np.array([int(np.where(r_all_powers==i)[0]) for i in r_sorted_powers[0:5]])
                r_peak_powers = r_all_powers[r_peak_indices]
                r_peak_periods = r_all_periods[r_peak_indices]
            
            else: 
                r_peak_powers = r_ls[0] 
                r_peak_periods = r_ls[1]
                r_all_powers = r_ls[2]
                r_all_periods = r_ls[3]
                r_faps = r_ls[4]
            
            rQ_results = quas_per(r_mjd, r_mag, r_magerr, r_peak_periods[0], sig_factor=1.25)
            rQ = rQ_results[0]
            r_resid_mags = rQ_results[1]
            
            rM = codyM(r_mag)
            
            # Other values
            rNu = sokolovskyNu(r_mag, r_magerr)
            rMean = np.mean(r_mag)
        
            # Run analysis for green data

            g_data_file = g_data_temp.replace('+++',id)
            g_data_df = pd.read_csv(g_data_file)
            g_mjd = np.array(g_data_df['mjd'])
            g_mag = np.array(g_data_df['mag'])
            g_magerr = np.array(g_data_df['magerr'])

            g_ls = lomb_scargle_routine(g_mjd, g_mag, g_magerr, minf, maxf)
            
            g_period_flag = g_ls[5]
            
            if g_period_flag == True: 
                g_faps = g_ls[4]
                
                g_cleaned_results = return_cleaned(g_ls[3], g_ls[2])
            
                g_all_powers = g_cleaned_results[1]
                g_all_periods = g_cleaned_results[0]
                
                g_sorted_powers = np.sort(g_all_powers)[::-1]
                g_peak_indices = np.array([int(np.where(g_all_powers==i)[0]) for i in g_sorted_powers[0:5]])
                g_peak_powers = g_all_powers[g_peak_indices]
                g_peak_periods = g_all_periods[g_peak_indices]
            
            else: 
                g_peak_powers = g_ls[0] 
                g_peak_periods = g_ls[1]
                g_all_powers = g_ls[2]
                g_all_periods = g_ls[3]
                g_faps = g_ls[4]        
            
            # GET PERIOD
            
            if r_peak_powers[0] > max(r_faps):
                per = r_peak_periods[0]
                
            elif g_peak_powers[0] > max(g_faps): 
                per = g_peak_periods[0]
                
            else: 
                per = np.nan
                
            if per > 249: 
                per = np.nan
            
            # Report results to results file
            
            if report_file != 'None' or report_file != 'none':
                with open(report_file,'a') as f:
                    id_string = id + ',' + str(ra) + ',' + str(dec)
                    results_string = str(per)+','+str(rQ)+','+str(rM)+','+str(rNu)+','+str(rMean) 
                    f.write(id_string+','+results_string+'\n')
            
            # Create Plot
            
            if out_type == 'None' or out_type == 'none':
                plt.clf()
                plt.close()

            elif out_type == 'show' or out_type == 'Show':
                make_plot(r_mjd, r_mag, r_magerr, r_peak_periods[0], r_all_periods, r_all_powers, r_faps, rQ, rM)
                plt.show()
                plt.clf()
            
            else:
                make_plot(r_mjd, r_mag, r_magerr, r_peak_periods[0], r_all_periods, r_all_powers, r_faps, rQ, rM)
                q_str = str(rQ) 
                q_str_split = q_str.split('.')
                decimal_portion = q_str_split[1]
                decimal_portion = decimal_portion[0:3]
                q_str = q_str_split[0]+'.'+decimal_portion
                actual_path = os.path.join(out_type,(q_str+':ID:'+id+'.png'))
                plt.savefig(actual_path, dpi=300, bbox_inches='tight') 
                plt.clf()
                plt.close()
                
            bar.next()
            
    bar.finish()

out_type = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/test/plots'
report_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/test/test_results.csv'
FAPs = [0.1, 0.05, 0.01]
key = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/key.csv'
r_data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/light_curves/+++_r.csv'
g_data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/light_curves/+++_r.csv'

june8thQ(key, r_data_temp, g_data_temp, min_per=0.5, max_per=250, false_alarm_levels=FAPs, out_type=out_type, report_file=report_file)