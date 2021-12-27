#Plot the interstellar reddening line
#irs was previously defined as 3.81
# utilizing ard3.81 ir slope, the ids GDR1_2163145984982389632,  2MASS_J20582381+4353114, and FHK_146 match great
# the range of slope values that could be consistent
# with a normal extinction law in the g vs g-r color magnitude diagram
# is 3.52 to 4.10



def slopes_demo_plot(ids, lc_dir, cmd_dir, slopes, intercepts, plot_dir):
    # Import(s)
    import os 
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    import matplotlib.ticker as ticker
    
    # Action
    
    # Create plot
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.figsize'] = (10, 4)
    
    fig, axs = plt.subplots(2, 2)
    
    # Edit some shared deets
    for i in range(2):
        axs[0, i].set_ylabel('r (mag)', fontsize=8)
        axs[0, i].set_xlabel('Time (MJD)', fontsize=8)
        
        lc_df = pd.read_csv(os.path.join(lc_dir, (ids[i]+'_r.csv')))
        mjds = np.array(lc_df['mjd'])
        mags = np.array(lc_df['mag'])
        magerrs = np.array(lc_df['magerr'])
        
        axs[0, i].errorbar(mjds, mags, magerrs, barsabove=False, elinewidth=0.2, ecolor='black', fmt='none', alpha=0.4)
        axs[0, i].scatter(mjds, mags, color='#408ee0', s=7, edgecolors='black', linewidths=0.3)
        axs[0, i].invert_yaxis()
        axs[0, i].set_xlabel('Date (MJD)')
        axs[0, i].set_ylabel('r (mag)')
        axs[1, i].margins(0.1, 0.1)
        axs[0, i].xaxis.set_major_locator(ticker.MultipleLocator(200))
        
        
        cmd_df = pd.read_csv(os.path.join(cmd_dir, (ids[i]+'_cmd.csv')))
        
        grs = np.array(cmd_df['color'])
        gs = np.array(cmd_df['g_mag'])
        med_g_err = np.median(cmd_df['g_magerr'])
        med_gr_err = np.median(cmd_df['color_err'])
        
        axs[1, i].margins(0.05, 0.1)
        
        axs[1, i].scatter(grs, gs, color='#408ee0', s=7, edgecolor='black', linewidths=0.3)        
        
        
        med_errs_x = axs[1, i].get_xlim()[0]
        med_errs_y = axs[1, i].get_ylim()[1]
        
        ## Best fit line
        slope = slopes[i]
        intercept = intercepts[i]
        
        x_bounds = np.array([(min(gs)-intercept), (max(gs)-intercept)])/slope
        
        bf_x = np.linspace(x_bounds[1],x_bounds[0],10)
        bf_y = (slope*bf_x)+intercept
        
        axs[1, i].plot(bf_x, bf_y, color='black', ls='--', label='Best Fit Slope', lw=0.5)
        
        axs[1, i].errorbar(x=med_errs_x, y=med_errs_y, xerr=med_gr_err, yerr=med_g_err, linewidth=0, elinewidth=1, capsize=1, label='Typical Errors', color='indianred')
        
        ####
        irs = 3.81 # interstellar redening slope
        max_y = max(gs)
        min_y = min(gs)
        med_y = float(np.median(gs))
        med_x = float(np.median(grs))

        max_x = ((max_y-med_y)/(irs))+med_x
        min_x = ((min_y-med_y)/(irs))+med_x

        x_vals = np.array([min_x,max_x,med_x])
        y_vals = (irs*(x_vals-med_x))+med_y

        axs[1, i].plot(x_vals,y_vals,color='indianred',label='Reddening Slope',ls='--', lw=0.5)
        ####
        
        axs[1, i].legend(loc='upper right', fontsize=6, fancybox=False, framealpha=0, edgecolor='black')
        axs[1, i].set_xlabel('g-r (mag)', fontsize=8)
        axs[1, i].set_ylabel('g (mag)', fontsize=8)
        axs[1, i].invert_yaxis()
        
        for j in range(2):
            axs[j, i].tick_params(axis='both', which='both', direction='in')
            axs[j, i].tick_params(which='both', bottom=True, top=True, left=True, right=True)
            axs[j, i].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
            axs[j, i].tick_params(axis='both', labelsize=8)
            axs[j, i].set_title(ids[i].replace('_', ' '), fontsize=8)
            axs[j, i].xaxis.set_minor_locator(AutoMinorLocator())
            axs[j, i].yaxis.set_minor_locator(AutoMinorLocator())
        
    plt.subplots_adjust(wspace=0.25, hspace=0.45)
    #plt.savefig(os.path.join(plot_dir, 'qpd_slopes_gallery_2.png'), bbox_inches='tight', dpi=600)
    plt.show()


### change all below after miles re-runs odr

### Bursters 
ids_b = ['FHK_32','FHK_142']
slopes_b = [1.9361835024995053, 0.985578160011938]
intercepts_b = [14.939790288192524, 14.667908951835173]
lc_dir = r'C:\Users\Research\Documents\GitHub\NAPYSOs\recovered\2.0\data\AUGUST_5th\light_curves'
cmd_dir = r'C:\Users\Research\Documents\GitHub\NAPYSOs\referee_changes\color-mag\reprocessing\cmd_files\cmd_files'
plot_dir = r''
slopes_demo_plot(ids=ids_b, lc_dir=lc_dir, cmd_dir=cmd_dir, slopes=slopes_b, intercepts=intercepts_b, plot_dir=plot_dir)
### dippers
ids_qpd = ['2MASS_J20582381+4353114','FHK_163']
slopes_qpd = [3.6602377592858373, 3.6589457101373295]
intercepts_qpd = [9.909792577788114, 8.465106990805705]
slopes_demo_plot(ids=ids_qpd, lc_dir=lc_dir, cmd_dir=cmd_dir, slopes=slopes_qpd, intercepts=intercepts_qpd, plot_dir=plot_dir)

