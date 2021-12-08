"""
This function is modified and should only be used for evaluating the window function things!!
"""

from matplotlib.pyplot import minorticks_on


def jumbled_q_m_and_ls(key,data_temp,min_per,max_per,rand_type,false_alarm_levels,out_type,report_file):
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
    import random
    from scipy.signal import find_peaks

    # Action

    key_df = pd.read_csv(key)
    ids = list(key_df['ID'])
    bar = Bar('Processing...',max=len(ids))

    if report_file != 'None' or report_file != 'none':
        with open(report_file,'a') as f:
            f.write('ID,POW1,PER1,POW2,PER2,POW3,PER3,POW4,PER4,90%_FAP,95%_FAP,99%_FAP'+'\n')

    '''
    We append values to the results file with each iteration rather than the usual 
    df to .csv file method with pandas so the majority of the work isn't lost if the 
    routine crashes prior to finishing. 
    '''   
    

    def create_plot():
        plt.rcParams['font.family'] = 'serif'
        fig, axs = plt.subplots(2,2)

        # Zoomed in periodogram plot
        axs[1, 0].plot(ls_periods,ls_powers,color='C0')
        axs[1, 0].set(xlim=(0.99, 1.01))      
        axs[1, 0].tick_params(axis='both',which='both',direction='in')
        axs[1, 0].tick_params(which='both',bottom=True,top=True,left=True,right=True)
        axs[1, 0].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        axs[1, 0].xaxis.set_minor_locator(AutoMinorLocator())
        axs[1, 0].yaxis.set_minor_locator(AutoMinorLocator()) 
        #axs[1, 0].set_xscale('log')
        axs[1, 0].set_xlabel('Period d')
        #axs[1, 0].set_yscale('log')
        axs[1, 0].set_ylabel('Power')
        #axs[1, 0].set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
        #axs[1, 0].set_ylim(0,1)
        axs[1, 0].set_title('Periodogram', fontsize=10)
        range_mask = np.logical_and(ls_periods>0.99, ls_periods<1.01)
        ranged_periods = ls_periods[range_mask]
        peak_locs, _ = find_peaks(ls_powers[range_mask], height=0.5)
        period_peaks = ranged_periods[peak_locs]
        
        for period in period_peaks: 
            axs[1, 0].axvline(x=period,label=str(round(period, 3))+' d',color='red')
                
        axs[1, 0].legend(loc='upper right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
        
        
        
        # Periodogram plot
        axs[1,1].plot(ls_periods,ls_powers,color='C0')
        axs[1,1].axvline(x=72,label='72 d',color='red')
        axs[1,1].axvline(x=91,label='91 d',color='red')
        axs[1,1].axvline(x=125,label='125 d',color='red')
        axs[1,1].axvline(x=224,label='224 d',color='red')
        #Plot FAP Levels
        colors = ['lightgrey','silver','darkgray','gray','dimgray']

        '''
        for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
            confidence_label = str(100*(1-false_alarm_levels[i]))[:-2]+'% FAP'
            axs[1,1].hlines(y=(faps[i]), xmin=min_per, xmax=max_per, color = colors[i],linestyles='--',
                label=confidence_label)
        '''

        axs[1,1].set_xscale('log')
        axs[1,1].set_xlabel('Period d')
        axs[1,1].set_xscale('log')
        axs[1,1].set_ylabel('Power')
        axs[1,1].set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
        axs[1,1].set_ylim(0,1)
        axs[1,1].set_xlim(min_per,max_per)
        axs[1,1].set_title('Periodogram', fontsize=10)
        axs[1,1].legend(loc='center',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
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

        plt.subplots_adjust(wspace=0.4, hspace=0.50)

    for id_index, id in enumerate(ids): 
        if id == 'GDR1_2162240708953684992': 
            # Get data 
            data_file = data_temp.replace('+++',id)
            data_df = pd.read_csv(data_file)
            mjd = np.array(data_df['mjd'])
            
            if rand_type==1: 
                mag = np.array(len(mjd)*[1])
            else: 
                mag = np.array(data_df['mag'])
                for i in range(10):
                    random.shuffle(mag)
                
                
                
            minf = 1/max_per
            maxf = 1/min_per

            # Periodogram of light curve
            periodogram = LombScargle(mjd,mag)
            ls_freqs,ls_powers = periodogram.autopower(method='fast',minimum_frequency=minf,maximum_frequency=maxf)
            faps = periodogram.false_alarm_level(false_alarm_levels)
            ls_periods = 1/ls_freqs
            
            # Find best period/power
            best_index = np.argmax(ls_powers)
            best_power = ls_powers[best_index]
            best_per = ls_periods[best_index]
            
            # find the rest
            peak_indices, _ = find_peaks(ls_powers,distance=9)
            peaks = np.sort(ls_powers[peak_indices])[::-1]
            second_index = np.where(ls_powers==peaks[1])
            second_per = float(ls_periods[second_index])
            second_pow = float(ls_powers[second_index])
            
            third_index = np.where(ls_powers==peaks[2])
            third_per = float(ls_periods[third_index])
            third_pow = float(ls_powers[third_index])
            
            fourth_index = np.where(ls_powers==peaks[3])
            fourth_per = float(ls_periods[fourth_index])
            fourth_pow = float(ls_powers[fourth_index])
            
    
            # Fold the light curve
            phased_dates = (np.mod(mjd, best_per))/best_per 
            phased_dates_cycle_2 = phased_dates + 1

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
                actual_path = out_type.replace('+++',id)
                plt.savefig(actual_path, dpi=300, format=out_type[-3:]) # Flexible save type (svg, png, etc.)
                plt.clf()
                plt.close()
            
            out_list = [best_power,best_per,second_pow,second_per,third_pow,third_per,fourth_pow,fourth_per]
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

# Ones for obs
key_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/key.csv'
data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
plot_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Jumble/ones_for_obs/plots/+++.png'
report_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Jumble/ones_for_obs/results.csv'
jumbled_q_m_and_ls(key=key_file,data_temp=data_temp,min_per=0.55,max_per=250,rand_type=1,false_alarm_levels=[0.1,0.05,0.01],out_type=plot_path,report_file=report_file)


# Jumbled obs 
'''
key_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/key.csv'
data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
plot_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Jumble/first_jumble/plots/+++.png'
report_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Jumble/first_jumble/results.csv'
jumbled_q_m_and_ls(key=key_file,data_temp=data_temp,min_per=0.55,max_per=250,rand_type=2,false_alarm_levels=[0.1,0.05,0.01],out_type=plot_path,report_file=report_file)
'''