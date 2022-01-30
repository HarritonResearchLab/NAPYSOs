def make_sup_plots(cmd_results, plot_dir):
    # Import(s)
    import os 
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator

    from tqdm import tqdm 
    
    # Action
    cmd_df = pd.read_csv(cmd_results)
    df = cmd_df
        
    new_ids = np.array(df['preferred_name'])
    ids = np.array(df['ID'])
    pers = np.array(df['period'])
    classes = np.array(df['primary_class'])
    qs = np.array(df['Q'])
    ms = np.array(df['M'])
    slopes = np.array(df['cmd_slope'])
    # use linear regression for intercepts 
    slope_angle_errors = np.array(df['angle_error'])
    
    '''
    ignore_ids = np.array([i.replace('.png', '') for i in os.listdir(plot_dir)])
    ids = np.setxor1d(ids, ignore_ids)
    '''
    for id, new_id, per, qm_class, q, m, slope, angle_err in tqdm(zip(
        ids, new_ids, pers, classes, qs, ms, slopes, slope_angle_errors)):
        
        
        lc_data = pd.read_csv('./recovered/2.0/data/AUGUST_5th/light_curves/'+id+'_r.csv') 
        cmd_data = pd.read_csv('./referee_changes/color-mag/reprocessing/cmd_files/cmd_files/'+id+'_cmd.csv') # fix
    
        
        mjds = np.array(lc_data['mjd'])
        mags = np.array(lc_data['mag'])
        magerrs = np.array(lc_data['magerr'])
        
        grs = np.array(cmd_data['color']) 
        gs = np.array(cmd_data['g_mag']) 
        med_g_err = np.median(cmd_data['g_magerr']) 
        med_gr_err = np.median(cmd_data['color_err']) 
        
        # Make plot
        plt.rcParams['font.family'] = 'serif'
        fig, axs = plt.subplots(1, 3, figsize=(12, 3))
        
        # Light curve
        ax = axs[0]
        
        title_str = new_id.replace('_',' ')+' '+ ' Q='+ str(round(q, 2)) 
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
        
        if pd.isnull(per) == False:
            phased_dates = (np.mod(mjds, per))/per 
            phased_dates_cycle_2 = phased_dates + 1
            ax.errorbar(phased_dates, mags, magerrs, barsabove=False, elinewidth=0.2, ecolor='black', fmt='none', alpha=0.4)
            ax.errorbar(phased_dates_cycle_2, mags, magerrs, barsabove=False, elinewidth=0.2, ecolor='black', fmt='none', alpha=0.4)
            ax.scatter(phased_dates, mags, color='#408ee0', s=7, edgecolors='black', linewidths=0.3)
            ax.scatter(phased_dates_cycle_2, mags, color='#408ee0', s=7, edgecolors='black', linewidths=0.3)
            ax.set_xlabel('Phase (P='+str(round(per,2))+' d)')
            ax.set_ylabel('r (mag)')
            ax.xaxis.set_major_locator(MultipleLocator(0.25))
            ax.invert_yaxis()
            
        
        # CMD data 
        ax = axs[2]
        
        ax.margins(0.05, 0.1)
        
        ax.scatter(grs, gs, color='#408ee0', s=7, edgecolor='black', linewidths=0.3)        
        
        med_errs_x = ax.get_xlim()[0]
        med_errs_y = ax.get_ylim()[1]
        
        ax.errorbar(x=med_errs_x, y=med_errs_y, xerr=med_gr_err, yerr=med_g_err, linewidth=0, elinewidth=1, capsize=1, label='Typical Errors', color='indianred')
        
        intercept = np.mean(gs - slope*grs)


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
        
        plot_path = plot_dir+new_id+'.png'
        plt.savefig(plot_path, dpi=200, bbox_inches='tight')
        plt.close(fig)

def find_missing(other_results, plot_dir):
    # Import(s)
    import os
    import numpy as np
    import pandas as pd 
    
    # Action
    ids = np.array(pd.read_csv(other_results)['ID'])
    ignore_ids = np.array([i.replace('.png', '') for i in os.listdir(plot_dir)])
    ids = np.setxor1d(ids, ignore_ids)
    print(ids)

def make_missing(ids, qs, ms, qm_classes, pers):
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator
    
    # Action
    for id, q, m, qm_class, per in zip(ids, qs, ms, qm_classes, pers):
        lc_data = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/'+id+'_r.csv')
        cmd_data = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/cmd/'+id+'.csv')
    
        mjds = np.array(lc_data['mjd'])
        mags = np.array(lc_data['mag'])
        magerrs = np.array(lc_data['magerr'])
        
        grs = np.array(cmd_data['g-r_val'])
        gs = np.array(cmd_data['g_mag'])
        med_g_err = np.median(cmd_data['g_magerr'])
        med_gr_err = np.median(cmd_data['g-r_magerr'])
        
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
        
        if pd.isnull(per) == False:
            phased_dates = (np.mod(mjds, per))/per 
            phased_dates_cycle_2 = phased_dates + 1
            ax.errorbar(phased_dates, mags, magerrs, barsabove=False, elinewidth=0.2, ecolor='black', fmt='none', alpha=0.4)
            ax.errorbar(phased_dates_cycle_2, mags, magerrs, barsabove=False, elinewidth=0.2, ecolor='black', fmt='none', alpha=0.4)
            ax.scatter(phased_dates, mags, color='#408ee0', s=7, edgecolors='black', linewidths=0.3)
            ax.scatter(phased_dates_cycle_2, mags, color='#408ee0', s=7, edgecolors='black', linewidths=0.3)
            ax.set_xlabel('Phase (P='+str(round(per,2))+' d)')
            ax.set_ylabel('r (mag)')
            ax.xaxis.set_major_locator(MultipleLocator(0.25))
            ax.invert_yaxis()
            
        
        # CMD data 
        ax = axs[2]
        
        ax.margins(0.05, 0.1)
        
        ax.scatter(grs, gs, color='#408ee0', s=7, edgecolor='black', linewidths=0.3)        
        
        med_errs_x = ax.get_xlim()[0]
        med_errs_y = ax.get_ylim()[1]
        
        ax.errorbar(x=med_errs_x, y=med_errs_y, xerr=med_gr_err, yerr=med_g_err, linewidth=0, elinewidth=1, capsize=1, label='Typical Errors', color='indianred')
        
        '''
        if 0 < angle_err < 10.0: 
        
            x_bounds = np.array([(min(gs)-intercept), (max(gs)-intercept)])/slope
            
            bf_x = np.linspace(x_bounds[1],x_bounds[0],10)
            bf_y = (slope*bf_x)+intercept
            
            ax.plot(bf_x, bf_y, color='black', ls='--', label='Best Fit Slope', lw=0.5)
        '''
        
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
        
        plot_path = plot_dir + '/'+id+'.png'
        plt.savefig(plot_path, dpi=200, bbox_inches='tight')
        plt.close(fig)

'''        
cmd_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Misc./supplimentary_figures/June_6th.csv'
other_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Misc./finalized_results/merged_results.csv'
plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Misc./supplimentary_figures/plots' 

make_sup_plots(cmd_results, other_results, plot_dir)
make_missing(['GDR1_2162123198648385024'], 
             [2.089131971415095],
             [-0.009262820775975348],
             ['u'],
             [27.117147485866408])
'''
 
cmd_results = './results/meta_data.csv'
plot_dir = './referee_changes/final_supp_plots/all_plots/'

make_sup_plots(cmd_results, plot_dir)