"""
This Function is the most up to date version of our Q/M/Lomb-Scargle routine 
"""

def q_m_and_ls(key,data_temp,min_per,max_per,pti,false_alarm_levels,out_type,report_file):
    ''' 
    ***
    Best working version of our lomb-scargle, flux-asymmetry, and quasi-periodicity code
    Certain sections were adopted and subsequently modified from code used in Bredall
    et al. (2020). Uses "hybrid" method rather than window function. 
    ***
    key: .csv file with a column of object ids labeled "ID"  
    data_temp: path to data file, with object id replaced with the string '+++' 
    min_per: minimum period 
    max_per: maximum period
    pti: i.e. Periods To Ignore, ranges of periods to ignore. E.g. 
            pti=['0.5:0.505','26:30'] ignores periods between 0.5-0.505
            and 26-30 days. 
    false_alarm_levels: FAPs, see documentation for astropy Lomb-Scargle
    out_type: 'Show' shows the plot, 'none' or 'None' skips creating the plot, 
            and anything else is used as the directory into which the plot
            should be saved. 
    '''
    # Import(s)
    from astropy.timeseries import LombScargle
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker
    #from progress.bar import Bar

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

        return quas_per, resid_mag

    # Action

    key_df = pd.read_csv(key)
    ids = list(key_df['ID'])
    #bar = Bar('Processing...',max=len(ids))

    if report_file != 'None' or report_file != 'none':
        with open(report_file,'a') as f:
            f.write('ID,POW1,PER1,POW2,PER2,POW3,PER3,POW4,PER4,OPOW1,OPER1,OPOW2,OPER2,OPOW3,OPER3,OPOW4,OPER4,OQ,Q,M,RPOW,RPER,90%_FAP,95%_FAP,99%_FAP'+'\n')

    '''
    We append values to the results file with each iteration rather than the usual 
    df to .csv file method with pandas so the majority of the work isn't lost if the 
    routine crashes prior to finishing. 
    '''   

    # Initializing these here to optimize speed
    harmonics = np.array(5*[1])/range(1,6)
    lower_harmonics = 0.98*harmonics
    upper_harmonics = 1.02*harmonics
    lph = 1/upper_harmonics # lph = Lower Period Harmonics (for lower bounds of harmonic period ranges)
    uph = 1/lower_harmonics # uph = Upper Period Harmonics (for upper bounds of harmonic period ranges)

    def per_is_good(period):
        harmonic = False
        for harmonic_idx in range(0,len(harmonics)):
            if lph[harmonic_idx] < period < uph[harmonic_idx]:
                harmonic = True
                break
        if harmonic==True: 
            return False
        else: 
            return True

    def create_plot():
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams.update({'font.size': 6})
        gs = gridspec.GridSpec(3,3)
        plt.suptitle('Q: '+str(round(q,3))+'; M: '+str(round(m,3)),fontsize=10)


        aliases = []

        if best_per > second_per:
            aliases.append(best_per/second_per)
        else:
            aliases.append(second_per/best_per)
        if best_per > third_per:
            aliases.append(best_per/third_per)
        else:
            aliases.append(third_per/best_per)
        if best_per > fourth_per:
            aliases.append(best_per/fourth_per)
        else:
            aliases.append(fourth_per/best_per)
        if second_per > third_per:
            aliases.append(second_per/third_per)
        else:
            aliases.append(third_per/second_per)
        if second_per > fourth_per:
            aliases.append(second_per/fourth_per)
        else:
            aliases.append(fourth_per/second_per)
        if third_per > fourth_per:
            aliases.append(third_per/fourth_per)
        else:
            aliases.append(fourth_per/third_per)

        x = 0
        alias_string = ""
        relationships = ["12", "13", "14", "23", "24", "34"]
        for a in aliases:
            if a >0.99 and a<1.01:
                alias_string = alias_string + (str(relationships[x]) + " ")
            if a >1.99 and a<2.01:
                alias_string = alias_string + (str(relationships[x]) + " ")
            if a >2.99 and a<3.01:
                alias_string = alias_string + (str(relationships[x]) + " ")
            if a >3.99 and a<4.01:
                alias_string = alias_string + (str(relationships[x]) + " ")
            if a >4.99 and a<5.01:
                alias_string = alias_string + (str(relationships[x]) + " ")
            x = x+1



        second_best = plt.subplot(gs[2,0])
        third_best = plt.subplot(gs[2, 1])
        fourth_best = plt.subplot(gs[2,2])
        periodogram = plt.subplot(gs[1, 0])
        best = plt.subplot(gs[1,1])
        resid = plt.subplot(gs[1, 2])
        lightcurve = plt.subplot(gs[0, 0:2])
        resid_periodogram = plt.subplot(gs[0, 2])
        plt.figtext(0.02,0.95,alias_string, fontdict={"fontsize": 20})

        #I need to create a lomb scargle of the residuals.
        #Also look for aliases in the results


        #ax[1,1]  is periodogram
        #ax[0,1] is best
        #ax[0,0] is light curve
        #ax[1, 0) is resid

        # Periodogram plot
        periodogram.plot(cleaned_periods,cleaned_powers,color='#408ee0')
        periodogram.plot(ls_period, ls_powers, color = "#40e093", lw = 0.3)



        for range_to_ignore in pti:
            range_to_ignore = range_to_ignore.split(':')
            lower_bound = float(range_to_ignore[0])
            upper_bound = float(range_to_ignore[1])
            periodogram.axvspan(xmin=lower_bound,xmax=upper_bound,alpha=1,color='white')
        
        #Plot FAP Levels
        colors = ['lightgrey','silver','darkgray','gray','dimgray']

        for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
            confidence_label = str(100*(1-false_alarm_levels[i]))[:-2]+'% FAP'
            periodogram.hlines(y=(faps[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
                label=confidence_label)
        
        def plot_ignored_ranges(): 
            for i in range(0,len(harmonics)):
                periodogram.axvspan(lph[i],uph[i],color='red',alpha=0.2)

            for i in pti:
                i_list = i.split(':')
                plt.axvspan(float(i_list[0]),float(i_list[1]),color='red',alpha=0.2)

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


        resid_periodogram.plot(cleaned_resid_periods, cleaned_resid_powers, color ="#408ee0")


        for range_to_ignore in pti:
            range_to_ignore = range_to_ignore.split(':')
            lower_bound = float(range_to_ignore[0])
            upper_bound = float(range_to_ignore[1])
            resid_periodogram.axvspan(xmin=lower_bound,xmax=upper_bound,alpha=1,color='white')

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


        # Folded light curve plot
        xlabel = 'Phase (P = '+str(round(best_per, 3))+' d)'
        best.set_xlabel(xlabel)
        best.set_ylabel('Mag')
        best.scatter(phased_dates, mag, s=2, color = "#408ee0")
        best.scatter(phased_dates_cycle_2, mag, s=2, color = "#408ee0")
        best.invert_yaxis()
        best.set_title("Best Folded Light Curve", fontsize=6)
        best.locator_params(axis='y', nbins=5)
        best.tick_params(axis='both',which='both',direction='in')
        best.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        best.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        best.xaxis.set_minor_locator(AutoMinorLocator())
        best.yaxis.set_minor_locator(AutoMinorLocator())


        #Second best folded light curve

        xlabel = 'Phase (P = '+str(round(second_per, 3))+' d)'
        second_best.set_xlabel(xlabel)
        second_best.set_ylabel('Mag')
        second_best.scatter(second_phased_dates, mag, s=2, color = "#408ee0")
        second_best.scatter(second_phased_dates_cycle_2, mag, s=2, color = "#408ee0")
        second_best.invert_yaxis()
        second_best.set_title("2nd best Folded Light Curve", fontsize=6)
        second_best.locator_params(axis='y', nbins=5)
        second_best.tick_params(axis='both',which='both',direction='in')
        second_best.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        second_best.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        second_best.xaxis.set_minor_locator(AutoMinorLocator())
        second_best.yaxis.set_minor_locator(AutoMinorLocator())

        xlabel = 'Phase (P = '+str(round(third_per, 3))+' d)'
        third_best.set_xlabel(xlabel)
        third_best.set_ylabel('Mag')
        third_best.scatter(third_phased_dates, mag, s=2, color = "#408ee0")
        third_best.scatter(third_phased_dates_cycle_2, mag, s=2, color = "#408ee0")
        third_best.invert_yaxis()
        third_best.set_title("3rd best Folded Light Curve", fontsize=6)
        third_best.locator_params(axis='y', nbins=5)
        third_best.tick_params(axis='both',which='both',direction='in')
        third_best.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        third_best.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        third_best.xaxis.set_minor_locator(AutoMinorLocator())
        third_best.yaxis.set_minor_locator(AutoMinorLocator())

        xlabel = 'Phase (P = '+str(round(fourth_per, 3))+' d)'
        fourth_best.set_xlabel(xlabel)
        fourth_best.set_ylabel('Mag')
        fourth_best.scatter(fourth_phased_dates, mag, s=2, color = "#408ee0")
        fourth_best.scatter(fourth_phased_dates_cycle_2, mag, s=2, color = "#408ee0")
        fourth_best.invert_yaxis()
        fourth_best.set_title("4th best Folded Light Curve", fontsize=6)
        fourth_best.locator_params(axis='y', nbins=5)
        fourth_best.tick_params(axis='both',which='both',direction='in')
        fourth_best.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        fourth_best.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        fourth_best.xaxis.set_minor_locator(AutoMinorLocator())
        fourth_best.yaxis.set_minor_locator(AutoMinorLocator())

        # Unfolded light curve plot
        lightcurve.errorbar(mjd, mag, yerr=magerr, lw=0,elinewidth=0.5, color = "#408ee0")
        lightcurve.scatter(mjd,mag,s=2, color = "#408ee0")
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
        resid.scatter(mjd,residuals,s=2, color = "#408ee0")
        resid.axhline(y=0,xmin=0,xmax=1,color='black',lw=0.75)
        resid.set_title('Residual Plot')
        resid.set_xlabel('MJD')
        resid.set_ylabel('Residual Mag',fontsize=6)
        resid.tick_params(axis='both',which='both',direction='in')
        resid.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        resid.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        resid.xaxis.set_minor_locator(AutoMinorLocator())
        resid.yaxis.set_minor_locator(AutoMinorLocator())

        plt.subplots_adjust(wspace=0.4, hspace=0.50)

    for id_index, id in enumerate(ids):

        if id_index > -1:
            # Get data 
            data_file = data_temp.replace('+++',id)
            data_df = pd.read_csv(data_file)
            mjd = np.array(data_df['mjd'])
            mag = np.array(data_df['mag'])
            magerr = np.array(data_df['magerr'])
            
            minf = 1/max_per
            maxf = 1/min_per

            # Periodogram of light curve
            periodogram = LombScargle(mjd,mag,magerr)
            ls_freqs,ls_powers = periodogram.autopower(method='fast',minimum_frequency=minf,maximum_frequency=maxf)
            faps = periodogram.false_alarm_level(false_alarm_levels)

            cleaned_powers = ls_powers
            cleaned_freqs = ls_freqs
            
        
            # CLEANING
            for f in cleaned_freqs:
                per = 1.0/f
                half_per = 0.5*per
                double_per = 2.0*per
                f_index = np.where(ls_freqs==f)

                if per_is_good(per)==True:
                    if per_is_good(half_per)==False or per_is_good(double_per)==False:
                        cleaned_freqs = np.delete(cleaned_freqs,f_index)
                        cleaned_powers = np.delete(cleaned_powers,f_index)
                else: 
                    cleaned_freqs = np.delete(cleaned_freqs,f_index)
                    cleaned_powers = np.delete(cleaned_powers,f_index)


            
            ### NEW METHOD FOR MASKING OUT OTHER RANGES
            cleaned_periods = 1/cleaned_freqs # 'cleaned' as in cleaned of significant window function periods
            ls_period = 1/ls_freqs
            for range_to_ignore in pti:
                range_to_ignore = range_to_ignore.split(':')
                lower_bound = float(range_to_ignore[0])
                upper_bound = float(range_to_ignore[1])
                indices_to_delete = np.logical_and(cleaned_periods<upper_bound,cleaned_periods>lower_bound)
                cleaned_periods = np.delete(arr=cleaned_periods,obj=indices_to_delete)
                cleaned_powers = np.delete(arr=cleaned_powers,obj=indices_to_delete)
            
            # Find best period/power
            best_index = np.argmax(cleaned_powers)
            best_power = cleaned_powers[best_index]
            best_per = cleaned_periods[best_index]

            ls_index = np.argmax(ls_powers)
            ls_power = ls_powers[ls_index]
            ls_per = ls_period[ls_index]


            # Find second period and its power
            from scipy.signal import find_peaks
            peak_indices, _ = find_peaks(cleaned_powers,distance=5)
            peaks = np.sort(cleaned_powers[peak_indices])[::-1]
            
            second_index = np.where(cleaned_powers==peaks[1])
            second_per = cleaned_periods[second_index]
            second_per = second_per[0]
            second_power = cleaned_powers[second_index]

            third_index = np.where(cleaned_powers==peaks[2])
            third_per = cleaned_periods[third_index]
            third_per = third_per[0]
            third_power = cleaned_powers[third_index]

            fourth_index = np.where(cleaned_powers==peaks[3])
            fourth_per = cleaned_periods[fourth_index]
            fourth_per = fourth_per[0]
            fourth_power = cleaned_powers[fourth_index]

            ls_peak_indices, _ = find_peaks(ls_powers,distance=5)
            ls_peaks = np.sort(ls_powers[ls_peak_indices])[::-1]

            ls_second_index = np.where(ls_powers==ls_peaks[1])
            ls_second_per = ls_period[ls_second_index]
            ls_second_per = ls_second_per[0]
            ls_second_power = ls_powers[ls_second_index]

            ls_third_index = np.where(ls_powers==ls_peaks[2])
            ls_third_per = ls_period[ls_third_index]
            ls_third_per = ls_third_per[0]
            ls_third_power = ls_powers[ls_third_index]

            ls_fourth_index = np.where(ls_powers==ls_peaks[3])
            ls_fourth_per = ls_period[ls_fourth_index]
            ls_fourth_per = ls_fourth_per[0]
            ls_fourth_power = ls_powers[ls_fourth_index]


            # Fold the light curve
            phased_dates = (np.mod(mjd, best_per))/best_per 
            phased_dates_cycle_2 = phased_dates + 1

            second_phased_dates = (np.mod(mjd, second_per))/second_per
            second_phased_dates_cycle_2 = second_phased_dates + 1

            third_phased_dates = (np.mod(mjd, third_per))/third_per
            third_phased_dates_cycle_2 = third_phased_dates + 1

            fourth_phased_dates = (np.mod(mjd, fourth_per))/fourth_per
            fourth_phased_dates_cycle_2 = fourth_phased_dates + 1


            # Calculate Q and get residuals plot
            qp_results = quas_per(mjd=mjd,mag=mag,magerr=magerr,per=best_per,sig_factor=1.25)
            q = qp_results[0]
            residuals = qp_results[1]

            ls_q_results = quas_per(mjd=mjd,mag=mag,magerr=magerr,per=ls_per,sig_factor=1.25)
            ls_q = ls_q_results[0]



            resid_periodogram = LombScargle(mjd, residuals, magerr)
            resid_freqs, resid_powers = resid_periodogram.autopower(method='fast',minimum_frequency=minf,maximum_frequency=maxf)
            resid_faps = resid_periodogram.false_alarm_level(false_alarm_levels)

            for f in resid_freqs:
                per = 1.0/f
                half_per = 0.5*per
                double_per = 2.0*per
                f_index = np.where(resid_freqs==f)

                if per_is_good(per)==True:
                    if per_is_good(half_per)==False or per_is_good(double_per)==False:
                        resid_freqs = np.delete(resid_freqs,f_index)
                        resid_powers = np.delete(resid_powers,f_index)
                else:
                    resid_freqs = np.delete(resid_freqs,f_index)
                    resid_powers = np.delete(resid_powers,f_index)

            cleaned_resid_powers = resid_powers
            cleaned_resid_freqs = resid_freqs

            cleaned_resid_periods = 1/cleaned_resid_freqs # 'cleaned' as in cleaned of significant window function periods

            for range_to_ignore in pti:
                range_to_ignore = range_to_ignore.split(':')
                lower_bound = float(range_to_ignore[0])
                upper_bound = float(range_to_ignore[1])
                indices_to_delete = np.logical_and(cleaned_resid_periods<upper_bound,cleaned_resid_periods>lower_bound)
                cleaned_resid_periods = np.delete(arr=cleaned_resid_periods,obj=indices_to_delete)
                cleaned_resid_powers = np.delete(arr=cleaned_resid_powers,obj=indices_to_delete)


            resid_index = np.argmax(cleaned_resid_powers)
            resid_power = cleaned_resid_powers[resid_index]
            resid_per = cleaned_resid_periods[resid_index]

            # Calculate m
            m = codyM(x=mag)

            # Plot
            if out_type == 'None' or out_type == 'none':
                plt.clf()
                plt.close()

            elif out_type == 'show' or out_type == 'Show':
                create_plot()
                plt.show()
                plt.clf()
            
            else:
                create_plot()
                
                q_str = str(q) 
                q_str_split = q_str.split('.')
                decimal_portion = q_str_split[1]
                decimal_portion = decimal_portion[0:3]
                q_str = q_str_split[0]+'.'+decimal_portion
                actual_path = out_type+'/'+q_str+':ID:'+id+'.png'
                plt.savefig(actual_path, dpi=300) 
                plt.clf()
                plt.close()
            
            out_list = [best_power, best_per,float(second_power),float(second_per), float(third_power), float(third_per), float(fourth_power), float(fourth_per), float(ls_power), float(ls_per),
                        float(ls_second_power), float(ls_second_per), float(ls_third_power), float(ls_third_per), float(ls_fourth_power), float(ls_fourth_per), ls_q,
                        q,m, float(resid_power), float(resid_per) ]
            
            for fap in faps: 
                out_list.append(fap)

            
            
            # Report results to result file
            if report_file != 'None' or report_file != 'none':
                with open(report_file,'a') as f:
                    result_string = id+','+','.join([str(val) for val in out_list])
                    f.write(result_string+'\n')




### RUN IT 

key_file = '/Users/s014605/Downloads/APRIL_7th 2/key.csv'
data_temp = '/Users/s014605/Downloads/APRIL_7th 2/light_curves/+++_g.csv'
plot_path = '/Users/s014605/Documents/HRL/Periodogram Plots_g/***:ID:+++.png'
report_file = '/Users/s014605/Documents/HRL/Periodogram Plots_g/results.csv'


q_m_and_ls(key=key_file,data_temp=data_temp,min_per=0.55,max_per=250,pti=['0.98:1.02','26:30'],false_alarm_levels=[0.1,0.05,0.01],out_type=plot_path,report_file=report_file)