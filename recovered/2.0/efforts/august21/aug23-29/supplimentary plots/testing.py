def test_plot(results_file, lc_dir, cmd_dir, plot_dir):
    # Import(s)
    import os
    import time
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator
    
    # Action 
    df = pd.read_csv(results_file)
    
    ids = np.array(df['ID'])
    slopes = np.array(df['odr slope'])
    intercepts = np.array(df['intercept'])
    angle_errs = np.array(df['angle error'])
    qs = np.array(df['Q'])
    ms = np.array(df['M'])
    periods = np.array(df['PER'])
    classes = np.array(df['primary_class'])
    
    for id in ids: 
        if id=='FHK_176': 
            id_loc = np.where(ids==id)
            slope = slopes[id_loc].item()
            intercept = intercepts[id_loc].item()
            angle_err = angle_errs[id_loc].item()
            q = qs[id_loc].item()
            m = ms[id_loc].item()
            per = periods[id_loc].item()
            qm_class = classes[id_loc].item()
                        
            lc_data = pd.read_csv(os.path.join(lc_dir, id+'_r.csv'))
            mjds = np.array(lc_data['mjd'])
            mags = np.array(lc_data['mag'])
            magerrs = np.array(lc_data['magerr'])
            
            cmd_data = pd.read_csv(os.path.join(cmd_dir, id+'.csv'))
            grs = np.array(cmd_data['g-r_val'])
            gs = np.array(cmd_data['g_mag'])
            med_g_err = np.median(cmd_data['g_magerr'])
            med_gr_err = np.median(cmd_data['g-r_magerr'])
            
            # Make plot
            # Make plot
            plt.rcParams['font.family'] = 'serif'
            fig, axs = plt.subplots(1, 3, figsize=(12, 3))
            
            # Light curve
            ax = axs[0]
            
            title_str = id.replace('_',' ')+' '+ ' Q='+ str(round(q, 2)) 
            title_str = title_str +'; M='+str(round(m, 2)) 
            title_str = title_str + ' ['+str(qm_class.upper())+']'
            
            fig.suptitle(title_str)
            
            ax.errorbar(mjds, mags, magerrs, barsabove=False, elinewidth=0.2, ecolor='black', fmt='none', alpha=0.4)
            ax.scatter(mjds, mags, color='#408ee0', s=7, edgecolors='black', linewidths=0.3)
            ax.invert_yaxis()
            ax.set_xlabel('Date (MJD)')
            ax.set_ylabel('r (mag)')
            
            
            # Folded light curve
            ax = axs[1]
            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            if pd.isnull(per) == False:
                phased_dates = (np.mod(mjds, per))/per 
                phased_dates_cycle_2 = phased_dates + 1
                ax.errorbar(phased_dates, mags, magerrs, barsabove=False, elinewidth=0.2, ecolor='black', fmt='none', alpha=0.4)
                ax.errorbar(phased_dates_cycle_2, mags, magerrs, barsabove=False, elinewidth=0.2, ecolor='black', fmt='none', alpha=0.4)
                ax.scatter(phased_dates, mags, color='#408ee0', s=7, edgecolors='black', linewidths=0.3)
                ax.scatter(phased_dates_cycle_2, mags, color='#408ee0', s=7, edgecolors='black', linewidths=0.3)
                ax.set_xlabel('Phase (P='+str(round(per,2))+' d)')
                ax.set_ylabel('r (mag)')
                ax.invert_yaxis()
                
            
            # CMD data 
            ax = axs[2]
            
            ax.margins(0.05, 0.1)
            
            ax.scatter(grs, gs, color='#408ee0', s=7, edgecolor='black', linewidths=0.3)        
            
            med_errs_x = ax.get_xlim()[0]
            med_errs_y = ax.get_ylim()[1]
            
            ax.errorbar(x=med_errs_x, y=med_errs_y, xerr=med_gr_err, yerr=med_g_err, linewidth=0, elinewidth=1, capsize=1, label='Typical Errors', color='indianred')
            
            if 0 < angle_err < 10.0: 
            
                x_bounds = np.array([(min(gs)-intercept), (max(gs)-intercept)])/slope
                
                bf_x = np.linspace(x_bounds[1],x_bounds[0],10)
                bf_y = (slope*bf_x)+intercept
                
                ax.plot(bf_x, bf_y, color='black', ls='--', label='Best Fit Slope', lw=0.5)
            
            ax.legend(loc='upper right', fontsize=6, fancybox=False, framealpha=0, edgecolor='black')
            ax.set_xlabel('g-r (mag)')
            ax.set_ylabel('g (mag)')
            ax.invert_yaxis()
            
            for ax in [axs[0], axs[1], axs[2]]:
                ax.tick_params(axis='both', which='both', direction='in')
                ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
                ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
                ax.tick_params(axis='both', labelsize=8)
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                
            axs[0].xaxis.set_major_locator(MultipleLocator(200))
                
            plt.subplots_adjust(wspace=0.35)
            #plt.show()
            plot_path = plot_dir + '/'+id+'.png'
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
            time.sleep(0.5)

results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv'
#lc_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/light_curves'
#cmd_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/cmds'
plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/aug23-29/supplimentary plots/plots'
id = 'GDR1_2162221673657763584'

lc_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/FHK_176'
cmd_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/FHK_176'

test_plot(results_file, lc_dir, cmd_dir, plot_dir)
    