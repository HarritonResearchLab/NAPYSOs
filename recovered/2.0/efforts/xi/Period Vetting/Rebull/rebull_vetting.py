def comparison_plots(our_results, rebull_results): 
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    
    # Action
    our_df = pd.read_csv(our_results)
    our_df = our_df[our_df['PER'].notna()]
    rebull_df = pd.read_csv(rebull_results)
    rebull_df = rebull_df[rebull_df[' period'].notna()]
    merged_df = pd.merge(left=our_df, right=rebull_df, sort=True, on='ID')
    ids = np.array(merged_df['ID'])
    
    our_pers = np.array(merged_df['PER'])
    rebull_pers = np.array(merged_df[' period'])
    
    flagged_ids = np.array([])
    our_flagged = np.array([])
    rebull_flagged = np.array([])
    
    for id, our_per, rebull_per in zip(ids, our_pers, rebull_pers):
        if np.abs((our_per-rebull_per)) > 0.1*(our_per): 
            our_flagged = np.append(our_flagged, our_per)
            rebull_flagged = np.append(rebull_flagged, rebull_per)
            flagged_ids = np.append(flagged_ids, id)
    
    def simple_comparison_plot(): 
        plt.rcParams['font.family']='serif'
        plt.rcParams['figure.figsize']=(9.5, 3.5)
        fig, axs = plt.subplots(1, 2)
        
        # All of ours versus all of hers
        axs[0].scatter(our_pers, rebull_pers, color='#408ee0', edgecolors='black', linewidths=0.5)
        axs[0].set_title('All Periodic Objects', fontsize='medium')
        axs[0].set_xscale('log')
        axs[0].set_yscale('log')
                
        axs[0].scatter(our_flagged, rebull_flagged, color='indianred', marker='o', edgecolors='black', linewidths=0.5, label='> 10% Difference')
        axs[0].legend(loc='lower right', fontsize='x-small', fancybox=False, edgecolor='black')
        
        
        for ax in [axs[0], axs[1]]:
            ax.tick_params(axis='both',which='both',direction='in')
            ax.tick_params(which='both',bottom=True,top=True,left=True,right=True)
            ax.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
            ax.set_xlabel('Our Period')
            ax.set_ylabel('Rebull Period')
                
        periodic_indices = np.where(merged_df['primary_class']=='p')
        our_p_pers = np.array(merged_df['PER'])[periodic_indices]
        rebull_p_pers = np.array(merged_df[' period'])[periodic_indices]
        
        axs[1].scatter(our_p_pers, rebull_p_pers, color='#408ee0', edgecolors='black', linewidths=0.5)
        axs[1].set_title('"P" Objects Only', fontsize='medium')
        axs[1].xaxis.set_minor_locator(AutoMinorLocator())    
        axs[1].yaxis.set_minor_locator(AutoMinorLocator())
        
        plt.subplots_adjust(wspace=0.4)
        plt.show()
    
    #simple_comparison_plot()
    
    def investigate_outliers(ids, data_temp, our_pers, rebull_pers, results_file):
        # Import(s)
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        # Definitions
        
        def cody_Q(mjd, mag, magerr, timescale, sig_factor):
            # Import(s)
            from astropy.convolution import Box1DKernel, convolve
            
            # Convert to arrays
            mjd = np.array(mjd)
            mag = np.array(mag)
            magerr = np.array(magerr)
                
            # Calculate sig
            sig = sig_factor*np.mean(magerr)

            # Create the residual curve
            phase = mjd % timescale
            mag = mag[np.argsort(phase)]

            # We use three periods and extract the middle to prevent edge effects
            three_periods = np.concatenate((mag, mag, mag))
            boxcar = Box1DKernel(len(mag) // 4)
            smooth_mag = convolve(three_periods, boxcar)
            smooth_mag = smooth_mag[np.size(mag):2*np.size(mag)]

            resid_mag = mag - smooth_mag

            q = float((np.nanstd(resid_mag) ** 2 - sig ** 2)/(np.nanstd(mag) ** 2 - sig ** 2))

            return q, resid_mag
        
        def make_diagnostic_plot(id, mjd, mag, our_per, rebull_per, our_q, rebull_q, plot_dir):
            # Import(s)
            import os 
            
            # Action
            plt.rcParams['font.family']='serif'
            plt.rcParams['figure.figsize']=(7, 4)
            fig, axs = plt.subplots(1, 2)
            
            for i in range(2):
                axs[i].tick_params(axis='both',which='both',direction='in')
                axs[i].tick_params(which='both',bottom=True,top=True,left=True,right=True)
                axs[i].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
            
            # Our period 
            phased_dates = (np.mod(mjd, our_per))/our_per 
            phased_dates_cycle_2 = phased_dates + 1
            axs[0].scatter(phased_dates, mag, color='#408ee0', s=3)
            axs[0].scatter(phased_dates_cycle_2, mag, color='#408ee0', s=3)
            axs[0].set_xlabel('Our Folded LC (P='+str(round(our_per, 2))+')')
            axs[0].set_ylabel('r (mag)')
            axs[0].set_title('Q: '+str(round(our_q, 2)), fontsize='medium')
            
            # Rebull period 
            phased_dates = (np.mod(mjd, rebull_per))/rebull_per 
            phased_dates_cycle_2 = phased_dates + 1
            axs[1].scatter(phased_dates, mag, color='#408ee0', s=3)
            axs[1].scatter(phased_dates_cycle_2, mag, color='#408ee0', s=3)
            axs[1].set_xlabel('Rebull Folded LC (P='+str(round(rebull_per, 2))+')')
            axs[1].set_ylabel('r (mag)')
            axs[1].set_title('Q: '+str(round(rebull_q, 2)), fontsize='medium')
            
            plt.subplots_adjust(wspace=0.4, hspace=0.3)
            
            plt.savefig(os.path.join(plot_dir, (id+'.png')))
            plt.clf()
            plt.close()
        
        def q_comparison_plot(our_qs, rebull_qs):
            plt.rcParams['font.family']='serif'
            plt.rcParams['figure.figsize']=(7, 3)
            
            plt.scatter(our_qs, rebull_qs, color='#408ee0', edgecolors='black', linewidths=0.5)
            plt.axhline(y=0.47, xmin=0, xmax=1, color='black')
            line_xs = np.array([min(our_qs), max(our_qs)])
            plt.plot(line_xs, (line_xs), ls='--', color='black')
            plt.xlabel('Our Q')
            plt.ylabel('Rebull Q')
            
            plt.tick_params(axis='both',which='both',direction='in')
            plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
            plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
            plt.show()
        
        
        our_qs = np.array([])
        rebull_qs = np.array([])
            
        for id, our_per, rebull_per in zip(ids, our_pers, rebull_pers):
            data_df = pd.read_csv(data_temp.replace('+++', id))
            mjds = np.array(data_df['mjd'])
            mags = np.array(data_df['mag'])
            magerrs = np.array(data_df['magerr'])
        
            our_q = cody_Q(mjds, mags, magerrs, our_per, 1.25)[0]
            rebull_q = cody_Q(mjds, mags, magerrs, rebull_per, 1.25)[0]
            
            our_qs = np.append(our_qs, our_q)
            rebull_qs = np.append(rebull_qs, rebull_q)
            
            plots_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Period Vetting/Rebull/investigating_discreps/plots'
            make_diagnostic_plot(id, mjds, mags, our_per, rebull_per, our_q, rebull_q, plots_dir)
        
        results_df = pd.DataFrame(list(zip(ids, our_pers, rebull_pers, our_qs, rebull_qs)), columns=['ID','OUR_PER','REB_PER', 'OUR_Q','REB_Q'])
        results_df.to_csv(results_file, index=False)
        
        q_comparison_plot(our_qs, rebull_qs)
    
    data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
    results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Period Vetting/Rebull/outliers.csv'
    investigate_outliers(flagged_ids, data_temp, our_flagged, rebull_flagged, results_file)
        
our_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Misc./finalized_results/merged_results.csv'
rebull_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Period Vetting/Rebull/rebull_period_assessment.csv'
comparison_plots(our_results, rebull_results)

def investigate_discreps(our_results, rebull_results): 
    # Import(s)
    import numpy as np
    import pandas as pd
    
    
