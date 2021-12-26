"""
This function is the q and ls funciton, but it uses the proposed Ansdell et al. 2018-Rodriguez et
al. 2017 "hybrid" method for accounting for aliases. This mfethod supersedes the theta hybrid_method.py
file. 
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
            and anything else is used as the path at which a copy of the plot 
            should be saved at. 
    '''
    # Import(s)
    from astropy.timeseries import LombScargle
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
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
        import pandas as pd
        from astropy.convolution import Box1DKernel, convolve

        # Action

        # Calculate sig (!!!)
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
    bar = Bar('Processing...',max=len(ids))

    if report_file != 'None' or report_file != 'none':
        with open(report_file,'a') as f:
            f.write('ID,POW1,PER1,POW2,PER2,Q,M,90%_FAP,95%_FAP,99%_FAP'+'\n')

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
        fig, axs = plt.subplots(2,2)
        fig.suptitle('Q: '+str(round(q,3))+'; M: '+str(round(m,3)),fontsize='medium')

        # Periodogram plot
        axs[1,1].plot(cleaned_periods,cleaned_powers,color='C0')

        for range_to_ignore in pti:
            range_to_ignore = range_to_ignore.split(':')
            lower_bound = float(range_to_ignore[0])
            upper_bound = float(range_to_ignore[1])
            axs[1,1].axvspan(xmin=lower_bound,xmax=upper_bound,alpha=1,color='white')
        
        #Plot FAP Levels
        colors = ['lightgrey','silver','darkgray','gray','dimgray']

        for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
            confidence_label = str(100*(1-false_alarm_levels[i]))[:-2]+'% FAP'
            axs[1,1].hlines(y=(faps[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
                label=confidence_label)
        
        def plot_ignored_ranges(): 
            for i in range(0,len(harmonics)):
                axs[1,1].axvspan(lph[i],uph[i],color='red',alpha=0.2)

            for i in pti:
                i_list = i.split(':')
                plt.axvspan(float(i_list[0]),float(i_list[1]),color='red',alpha=0.2)

        axs[1,1].set_xscale('log')
        axs[1,1].set_xlabel('Period d')
        axs[1,1].set_xscale('log')
        axs[1,1].set_ylabel('Power')
        axs[1,1].set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
        axs[1,1].set_ylim(0,1)
        axs[1,1].set_xlim(min_per,max_per)
        axs[1,1].set_title('Periodogram', fontsize=10)
        axs[1,1].legend(loc='upper right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
        axs[1,1].tick_params(axis='both',which='both',direction='in')
        axs[1,1].tick_params(which='both',bottom=True,top=True,left=True,right=True)
        axs[1,1].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())

        # Folded light curve plot
        xlabel = 'Phase (P = '+str(round(best_per, 3))+' d)'
        axs[0,1].set_xlabel(xlabel)
        axs[0,1].set_ylabel('Mag')
        axs[0,1].scatter(phased_dates, mag, s=2)
        axs[0,1].scatter(phased_dates_cycle_2, mag, s=2, c='C0')
        axs[0,1].invert_yaxis()
        axs[0,1].set_title("Folded Light Curve", fontsize=10)
        axs[0,1].locator_params(axis='y', nbins=5)
        axs[0,1].tick_params(axis='both',which='both',direction='in')
        axs[0,1].tick_params(which='both',bottom=True,top=True,left=True,right=True)
        axs[0,1].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
        axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())

        # Unfolded light curve plot
        axs[0,0].errorbar(mjd, mag, yerr=magerr, lw=0,elinewidth=0.5)
        axs[0,0].scatter(mjd,mag,s=2)
        axs[0,0].invert_yaxis()
        axs[0,0].set_ylabel('Mag')
        axs[0,0].set_xlabel('MJD')
        axs[0,0].set_title('Light Curve', fontsize=10)
        axs[0,0].locator_params(axis='y', nbins=5)
        axs[0,0].tick_params(axis='both',which='both',direction='in')
        axs[0,0].tick_params(which='both',bottom=True,top=True,left=True,right=True)
        axs[0,0].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
        axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())

        # Residuals plot
        axs[1,0].scatter(mjd,residuals,s=2)
        axs[1,0].axhline(y=0,xmin=0,xmax=1,color='black',lw=0.75)
        axs[1,0].set_title('Residual Plot')
        axs[1,0].set_xlabel('MJD')
        axs[1,0].set_ylabel('Residual Mag',fontsize=10)
        axs[1,0].tick_params(axis='both',which='both',direction='in')
        axs[1,0].tick_params(which='both',bottom=True,top=True,left=True,right=True)
        axs[1,0].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
        axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())

        plt.subplots_adjust(wspace=0.4, hspace=0.50)

    for id_index, id in enumerate(ids): 
        if id_index>244: 
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
            
        
            # CLEANING
            for f in ls_freqs: 
                per = 1.0/f
                half_per = 0.5*per
                double_per = 2.0*per
                f_index = np.where(ls_freqs==f)

                if per_is_good(per)==True:
                    if per_is_good(half_per)==False or per_is_good(double_per)==False:
                        ls_freqs = np.delete(ls_freqs,f_index)
                        ls_powers = np.delete(ls_powers,f_index)
                else: 
                    ls_freqs = np.delete(ls_freqs,f_index)
                    ls_powers = np.delete(ls_powers,f_index)
                        
            cleaned_powers = ls_powers
            cleaned_freqs = ls_freqs

            
            ### NEW METHOD FOR MASKING OUT OTHER RANGES
            cleaned_periods = 1/cleaned_freqs # 'cleaned' as in cleaned of significant window function periods
            
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

            # Find second period and its power
            from scipy.signal import find_peaks
            peak_indices, _ = find_peaks(cleaned_powers,distance=5)
            peaks = np.sort(cleaned_powers[peak_indices])[::-1]
            
            second_index = np.where(cleaned_powers==peaks[1])
            second_per = cleaned_periods[second_index]
            second_power = cleaned_powers[second_index]

            # Fold the light curve
            phased_dates = (np.mod(mjd, best_per))/best_per 
            phased_dates_cycle_2 = phased_dates + 1

            # Calculate Q and get residuals plot
            qp_results = quas_per(mjd=mjd,mag=mag,magerr=magerr,per=best_per,sig_factor=1.25)
            q = qp_results[0]
            residuals = qp_results[1]

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
                actual_path = out_type.replace('***',q_str)
                actual_path = actual_path.replace('+++',id)
                plt.savefig(actual_path, dpi=300, format=out_type[-3:]) # Flexible save type (svg, png, etc.)
                plt.clf()
                plt.close()
            
            out_list = [best_power, best_per,float(second_power),float(second_per),q,m]
            for fap in faps: 
                out_list.append(fap)

            
            
            # Report results to result file
            if report_file != 'None' or report_file != 'none':
                with open(report_file,'a') as f:
                    result_string = id+','+','.join([str(val) for val in out_list])
                    f.write(result_string+'\n')

            bar.next()
    
    bar.finish()

### RUN IT 

key_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/key.csv'
data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
plot_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/third/plots/***:ID:+++.png'
report_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/third/results.csv'

q_m_and_ls(key=key_file,data_temp=data_temp,min_per=0.55,max_per=250,pti=['0.95:1.05','26:30'],false_alarm_levels=[0.1,0.05,0.01],out_type=plot_path,report_file=report_file)