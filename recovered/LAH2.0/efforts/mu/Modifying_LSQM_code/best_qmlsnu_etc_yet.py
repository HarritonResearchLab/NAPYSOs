



def nearly_everything(key,r_data_temp, g_data_temp,min_per,max_per,false_alarm_levels,r_out_type,g_out_type,flag_dir, sig_g_dir, report_file):
    ''' 
    ***
    Best working version of our lomb-scargle, flux-asymmetry, and quasi-periodicity code
    Certain sections were adopted and subsequently modified from code used in Bredall
    et al. (2020). This function has been modified so it doesn't ignore or clean stuff; it just 
    flags if there is no dominant red period but a signifiant green one or if the peak red
    period is within the canidate aliase ranges. This function also does all the analysis for both
    the red and green data. 
    ***
    key: .csv file with a column of object ids labeled "ID"  
    r_data_temp: path to g data file, with object id replaced with the string '+++' 
    g_data_temp: path to g data file, "" 
    min_per: minimum period. Set to 0.5
    max_per: maximum period. Set to 250
    false_alarm_levels: FAPs, see documentation for astropy Lomb-Scargle.
                        Set to [0.1, 0.05, 0.01]
    r_out_type: 'Show' shows the plot, 'none' or 'None' skips creating the plot, 
            and anything else is used as the directory into which the plot
            should be saved. Note that this should be the directory, not the usual 
            file path with the id replaced with '+++'. 
    g_out_type: Same as for red, but for the green plots. 
    flag_dir: the directory into which the file containing the flagged objects should be
              deposited. 
    sig_g_dir: the directory into which the file containing the sig. g 
               periods should be deposited. 
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
        harmonics = np.array([0.5,1,2]) # What's the current range around 1? 
        lph = 0.98*harmonics
        uph = 1.02*harmonics
        
        def check_per(period):
            harmonic = False
            for harmonic_idx in range(0,len(harmonics)):
                if lph[harmonic_idx] < period < uph[harmonic_idx]:
                    harmonic = True
                    break
            if harmonic==True: 
                return False
            else: 
                return True
        
        # Periodogram of light curve
        periodogram = LombScargle(mjd,mag,magerr)
        all_freqs, all_powers = periodogram.autopower(method='fast',minimum_frequency=minf,maximum_frequency=maxf)
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
            half_per = 0.5*period
            double_per = 2.0*period
            for per in [period, half_per, double_per]:
                for middle in [0.5, 1, 2]:
                    if 0.98*middle < per < 1.02*middle:
                        return True 
                        break
                    else: 
                        if 26 < per < 30: 
                            return True
                            break
                        else: 
                            return True 
        
        return peak_powers, peak_periods, all_powers, all_periods, faps, flag(peak_periods[0]) 
     
    def diagnostic_plot(Q, M, mjds, mags, magerrs, resid_mags, peak_periods, all_periods, all_powers, faps):
        # Import(s)
        from matplotlib.ticker import AutoMinorLocator
        import matplotlib.gridspec as gridspec
        import matplotlib.ticker
        from astropy.timeseries import LombScargle
        
        def fold_lc(mjds, per):
            phased_dates = (np.mod(mjds, per))/per 
            phased_dates_cycle_2 = phased_dates + 1
            return phased_dates, phased_dates_cycle_2
        
        # Action
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams.update({'font.size': 6})
        gs = gridspec.GridSpec(3,3)
        plt.suptitle('Q: +'+str(round(Q,3))+'; M: '+str(round(M, 3)),fontsize=10)

        second_best = plt.subplot(gs[2,0])
        third_best = plt.subplot(gs[2, 1])
        fourth_best = plt.subplot(gs[2,2])
        periodogram = plt.subplot(gs[1, 0])
        best = plt.subplot(gs[1,1])
        resid = plt.subplot(gs[1, 2])
        lightcurve = plt.subplot(gs[0, 0:2])
        resid_periodogram = plt.subplot(gs[0, 2])


        # Periodogram plot
        periodogram.plot(all_periods,all_powers,color='#408ee0')
        
        #Plot FAP Levels
        colors = ['lightgrey','silver','darkgray','gray','dimgray']

        for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
            confidence_label = str(100*(1-false_alarm_levels[i]))[:-2]+'% FAP'
            periodogram.hlines(y=(faps[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
                label=confidence_label)

        periodogram.set_xscale('log')
        periodogram.set_xlabel('Period d')
        periodogram.set_xscale('log')
        periodogram.set_ylabel('Power')
        periodogram.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
        periodogram.set_ylim(0,1)
        periodogram.set_xlim(min_per,max_per)
        periodogram.set_title('Periodogram', fontsize=6)
        periodogram.tick_params(axis='both',which='both',direction='in')
        periodogram.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        periodogram.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(0.158, 0.2511, 0.398,0.63 ))
        periodogram.xaxis.set_minor_locator(locmaj)
        periodogram.yaxis.set_minor_locator(AutoMinorLocator())
        periodogram.legend(loc='upper right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)

        # Periodogram on residuals
        periodogram = LombScargle(mjds, resid_mags)
        resid_ls_freqs,resid_ls_powers = periodogram.autopower(method='fast',minimum_frequency=minf,maximum_frequency=maxf)
        resid_faps = periodogram.false_alarm_level(false_alarm_levels)
        resid_pers = 1/resid_ls_freqs
        
        resid_periodogram.plot(resid_pers, resid_ls_powers, color ="#408ee0")

        for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is
            confidence_label = str(100*(1-false_alarm_levels[i]))[:-2]+'% FAP'
            resid_periodogram.hlines(y=(resid_faps[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
                label=confidence_label)

        resid_periodogram.set_xscale('log')
        resid_periodogram.set_xlabel('Period d')
        resid_periodogram.set_xscale('log')
        resid_periodogram.set_ylabel('Power')
        resid_periodogram.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
        resid_periodogram.set_xlim(min_per,max_per)
        resid_periodogram.set_title('Residual Periodogram', fontsize=6)
        resid_periodogram.legend(loc='upper right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
        resid_periodogram.tick_params(axis='both',which='both',direction='in')
        resid_periodogram.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        resid_periodogram.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(0.158, 0.2511, 0.398,0.63 ))
        resid_periodogram.xaxis.set_minor_locator(locmaj)
        resid_periodogram.yaxis.set_minor_locator(AutoMinorLocator())

        # Best folded light curve plot 
        best_folded = fold_lc(mjds, peak_periods[0])
        phased_dates = best_folded[0]
        phased_dates_cycle_2 = best_folded[1]
        
        xlabel = 'Phase (P = '+str(round(peak_periods[0], 3))+' d)'
        best.set_xlabel(xlabel)
        best.set_ylabel('Mag')
        best.scatter(phased_dates, mags, s=2, color = "#408ee0")
        best.scatter(phased_dates_cycle_2, mags, s=2, color = "#408ee0")
        best.invert_yaxis()
        best.set_title("Best Folded LC", fontsize=6)
        best.locator_params(axis='y', nbins=5)
        best.tick_params(axis='both',which='both',direction='in')
        best.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        best.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        best.xaxis.set_minor_locator(AutoMinorLocator())
        best.yaxis.set_minor_locator(AutoMinorLocator())

        #Second best folded light curve
        second_best_folded = fold_lc(mjds, peak_periods[1])
        second_phased_dates = second_best_folded[0]
        second_phased_dates_cycle_2 = second_best_folded[1]
        xlabel = 'Phase (P = '+str(round(peak_periods[1], 3))+' d)'
        second_best.set_xlabel(xlabel)
        second_best.set_ylabel('Mag')
        second_best.scatter(second_phased_dates, mags, s=2, color = "#408ee0")
        second_best.scatter(second_phased_dates_cycle_2, mags, s=2, color = "#408ee0")
        second_best.invert_yaxis()
        second_best.set_title("2nd Best Folded LC", fontsize=6)
        second_best.locator_params(axis='y', nbins=5)
        second_best.tick_params(axis='both',which='both',direction='in')
        second_best.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        second_best.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        second_best.xaxis.set_minor_locator(AutoMinorLocator())
        second_best.yaxis.set_minor_locator(AutoMinorLocator())

        # Third Best Folded Lightcurve
        third_best_folded = fold_lc(mjds, peak_periods[2])
        third_phased_dates = third_best_folded[0]
        third_phased_dates_cycle_2 = third_best_folded[1]
        xlabel = 'Phase (P = '+str(round(peak_periods[2], 3))+' d)'
        third_best.set_xlabel(xlabel)
        third_best.set_ylabel('Mag')
        third_best.scatter(third_phased_dates, mags, s=2, color = "#408ee0")
        third_best.scatter(third_phased_dates_cycle_2, mags, s=2, color = "#408ee0")
        third_best.invert_yaxis()
        third_best.set_title("3rd Best Folded LC", fontsize=6)
        third_best.locator_params(axis='y', nbins=5)
        third_best.tick_params(axis='both',which='both',direction='in')
        third_best.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        third_best.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        third_best.xaxis.set_minor_locator(AutoMinorLocator())
        third_best.yaxis.set_minor_locator(AutoMinorLocator())

        # Fourth Best Folded Lightcurve 
        fourth_best_folded = fold_lc(mjds, peak_periods[3])
        fourth_phased_dates = fourth_best_folded[0]
        fourth_phased_dates_cycle_2 = fourth_best_folded[1]
        xlabel = 'Phase (P = '+str(round(peak_periods[3], 3))+' d)'
        fourth_best.set_xlabel(xlabel)
        fourth_best.set_ylabel('Mag')
        fourth_best.scatter(fourth_phased_dates, mags, s=2, color = "#408ee0")
        fourth_best.scatter(fourth_phased_dates_cycle_2, mags, s=2, color = "#408ee0")
        fourth_best.invert_yaxis()
        fourth_best.set_title("4th Best Folded LC", fontsize=6)
        fourth_best.locator_params(axis='y', nbins=5)
        fourth_best.tick_params(axis='both',which='both',direction='in')
        fourth_best.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        fourth_best.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        fourth_best.xaxis.set_minor_locator(AutoMinorLocator())
        fourth_best.yaxis.set_minor_locator(AutoMinorLocator())

        # Unfolded light curve plot
        lightcurve.errorbar(mjds, mags, yerr=magerrs, lw=0,elinewidth=0.5, color = "#408ee0")
        lightcurve.scatter(mjds,mags,s=2, color = "#408ee0")
        lightcurve.invert_yaxis()
        lightcurve.set_ylabel('Mag')
        lightcurve.set_xlabel('MJD')
        lightcurve.set_title('Light Curve', fontsize=6)
        lightcurve.locator_params(axis='y', nbins=5)
        lightcurve.tick_params(axis='both',which='both',direction='in')
        lightcurve.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        lightcurve.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        lightcurve.xaxis.set_minor_locator(AutoMinorLocator())
        lightcurve.yaxis.set_minor_locator(AutoMinorLocator())

        # Residuals plot
        resid.scatter(mjds,resid_mags,s=2, color = "#408ee0")
        resid.axhline(y=0,xmin=0,xmax=1,color='black',lw=0.75)
        resid.set_title('Residual Plot', fontsize=6)
        resid.set_xlabel('MJD')
        resid.set_ylabel('Residual Mag',fontsize=6)
        resid.tick_params(axis='both',which='both',direction='in')
        resid.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        resid.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        resid.xaxis.set_minor_locator(AutoMinorLocator())
        resid.yaxis.set_minor_locator(AutoMinorLocator())

        plt.subplots_adjust(wspace=0.4, hspace=0.50)
    
    # Action

    key_df = pd.read_csv(key)
    ids = list(key_df['ID'])
    ras = np.array(key_df['RA'])
    decs = np.array(key_df['DEC'])
    
    bar = Bar('Processing...', max=len(ids))

    if report_file != 'None' or report_file != 'none':
        with open(report_file,'a') as f:
            r_columns = 'rPOW1,rPER1,rPOW2,rPER2,rPOW3,rPER3,rPOW4,rPER4,rQ,rM,rNu,rMean,r90%_FAP,r95%_FAP,r99%_FAP,flag'
            g_columns = r_columns.replace('rM,rNu,rMean,','')
            g_columns = g_columns.replace('r','g')
            g_columns = g_columns.replace(',flag','')
            f.write('ID,RA,DEC,'+r_columns+','+g_columns+'\n')
           
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
            
            r_peak_powers = r_ls[0] 
            r_peak_periods = r_ls[1]
            r_all_powers = r_ls[2]
            r_all_periods = r_ls[3]
            r_faps = r_ls[4]
            period_flag = r_ls[5]
            
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
            
            g_peak_powers = g_ls[0] 
            g_peak_periods = g_ls[1]
            g_all_powers = g_ls[2]
            g_all_periods = g_ls[3]
            g_faps = g_ls[4]           
            
            gQ_results = quas_per(g_mjd, g_mag, g_magerr, g_peak_periods[0], sig_factor=1.25)
            gQ = gQ_results[0]
            g_resid_mags = gQ_results[1]
            
            secondFlag = False
            # Evaluate conditions for flagging 
            if period_flag == True and rQ < 0.5: 
                flag = True
            else: 
                flag = False
            if rQ>0.45 and gQ < 0.5:
                secondFlag = True
            else: 
                secondFlag = False 
        
            if report_file != 'None' or report_file != 'none':
                with open(report_file,'a') as f:
                    
                    if flag or secondFlag == True:
                        finalFlag = True
                    else: 
                        finalFlag = False

                    id_string = id + ',' + str(ra) + ',' + str(dec)
                    
                    r_results = ','.join([(str(r_peak_powers[i])+','+str(r_peak_periods[i])) for i in range(4)])+','
                    r_results = r_results + str(rQ) + ',' + str(rM) + ',' + str(rNu) + ',' + str(rMean) + ','
                    r_results = r_results + ','.join([str(fap) for fap in r_faps[0:4]]) + ',' + str(finalFlag)
                    
                    g_results = ','.join([(str(g_peak_powers[i])+','+str(g_peak_periods[i])) for i in range(4)])+','
                    g_results = g_results + str(gQ) + ',' + ','.join([str(fap) for fap in g_faps[0:4]])
                    
                    f.write(id_string + ',' + r_results + ',' + g_results + '\n')
            
            # Create Plots
            
            # r Plot
            if r_out_type == 'None' or r_out_type == 'none':
                plt.clf()
                plt.close()

            elif r_out_type == 'show' or r_out_type == 'Show':
                diagnostic_plot(rQ, rM, r_mjd, r_mag, r_magerr, r_resid_mags, r_peak_periods, r_all_periods, r_all_powers, r_faps)
                plt.show()
                plt.clf()
            
            else:
                diagnostic_plot(rQ, rM, r_mjd, r_mag, r_magerr, r_resid_mags, r_peak_periods, r_all_periods, r_all_powers, r_faps)
                q_str = str(rQ) 
                q_str_split = q_str.split('.')
                decimal_portion = q_str_split[1]
                decimal_portion = decimal_portion[0:3]
                q_str = q_str_split[0]+'.'+decimal_portion
                actual_path = os.path.join(r_out_type,(q_str+':ID:'+id+'.png'))
                flag_path = os.path.join(flag_dir,(q_str+':ID:'+id+'.png'))
                sig_g_path = os.path.join(sig_g_dir, (q_str+':ID:'+id+'.png'))
                plt.savefig(actual_path, dpi=300) 
                if flag == True: 
                    plt.savefig(flag_path, dpi=300)
                if secondFlag == True: 
                    plt.savefig(sig_g_path, dpi=300)
                plt.clf()
                plt.close()
                
            # g Plot
            if g_out_type == 'None' or g_out_type == 'none':
                plt.clf()
                plt.close()

            elif g_out_type == 'show' or g_out_type == 'Show':
                diagnostic_plot(gQ, 0, g_mjd, g_mag, g_magerr, g_resid_mags, g_peak_periods, g_all_periods, g_all_powers, g_faps)
                plt.show()
                plt.clf()
            
            else:
                diagnostic_plot(gQ, 0, g_mjd, g_mag, g_magerr, g_resid_mags, g_peak_periods, g_all_periods, g_all_powers, g_faps)
                q_str = str(gQ) 
                q_str_split = q_str.split('.')
                decimal_portion = q_str_split[1]
                decimal_portion = decimal_portion[0:3]
                q_str = q_str_split[0]+'.'+decimal_portion
                actual_path = os.path.join(g_out_type, (q_str+':ID:'+id+'.png'))
                plt.savefig(actual_path, dpi=300) 
                plt.clf()
                plt.close()
                
            bar.next()
            
    bar.finish()

r_out_type = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/qm_again/plots/r_plots'
g_out_type = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/qm_again/plots/g_plots'
flag_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/qm_again/plots/flagged'
report_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/qm_again/qm_plus_results.csv'
FAPs = [0.1, 0.05, 0.01]
key = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/key.csv'
r_data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/light_curves/+++_r.csv'
g_data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/light_curves/+++_r.csv'

nearly_everything(key, r_data_temp, g_data_temp, min_per=0.5, max_per=250, false_alarm_levels=FAPs, r_out_type=r_out_type, g_out_type=g_out_type, flag_dir = flag_dir, report_file=report_file)