def slopes_demo_plot(ids, lc_dir, cmd_dir, slopes, intercepts, plot_dir):
    # Import(s)
    import os 
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    
    # Action
    
    # Create plot
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.figsize'] = (11, 4)
    
    fig, axs = plt.subplots(2, 3)
    
    # Edit some shared deets
    for i in range(3):
        axs[0, i].set_ylabel('r (mag)', fontsize=8)
        axs[0, i].set_xlabel('Time (MJD)', fontsize=8)
        
        lc_df = pd.read_csv(os.path.join(lc_dir, (ids[i]+'_r.csv')))
        mjds = np.array(lc_df['mjd'])
        mags = np.array(lc_df['mag'])
        magerrs = np.array(lc_df['magerr'])
        
        axs[0, i].errorbar(mjds, mags, magerrs, color='#408ee0', barsabove=False, marker='o', ms=7, lw=0, elinewidth=0.2, ecolor='black')
        #axs[0, i].scatter(mjds, mags, color='#408ee0', s=7, edgecolor='black', linewidths=0.35)
        axs[0, i].invert_yaxis()
        axs[0, i].set_xlabel('Date (MJD)')
        axs[0, i].set_ylabel('r (mag)')
        axs[1, i].margins(0.1, 0.1)
        
        
        cmd_df = pd.read_csv(os.path.join(cmd_dir, (ids[i]+'.csv')))
        
        grs = np.array(cmd_df['g-r_val'])
        gs = np.array(cmd_df['g_mag'])
        med_g_err = np.median(cmd_df['g_magerr'])
        med_gr_err = np.median(cmd_df['g-r_magerr'])
        
        axs[1, i].margins(0.05, 0.1)
        
        axs[1, i].scatter(grs, gs, color='#408ee0', s=0.75)#, edgecolor='black', linewidths=0.1)        
        
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
        
    plt.subplots_adjust(wspace=0.35, hspace=0.5)
    #plt.show()
    plt.savefig(os.path.join(plot_dir, 'qpd_slopes_gallery.png'), bbox_inches='tight', dpi=600)

def slopes_demo_plot_2(ids, lc_dir, cmd_dir, slopes, intercepts, plot_dir):
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
        
        
        cmd_df = pd.read_csv(os.path.join(cmd_dir, (ids[i]+'.csv')))
        
        grs = np.array(cmd_df['g-r_val'])
        gs = np.array(cmd_df['g_mag'])
        med_g_err = np.median(cmd_df['g_magerr'])
        med_gr_err = np.median(cmd_df['g-r_magerr'])
        
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
    plt.savefig(os.path.join(plot_dir, 'qpd_slopes_gallery_2.png'), bbox_inches='tight', dpi=600)

def find_objects(lc_dir, odr_results, qm_results, plot_dir):
    # Import(s)
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Action
    #key_df = pd.read_csv(key)
    #key_ids = np.array(key_df['ID'])
    
    slopes_df = pd.read_csv(odr_results)
    slope_ids = np.array(slopes_df['ID'])
    slopes = np.array(slopes_df['SLOPE_ANGLE'])
    slope_errors = np.array(slopes_df['SLOPE_ANGLE_ERROR'])
    
    sig_slopes_indices = np.intersect1d(np.where(slope_errors<10.0),np.where(slopes>0))

    slope_ids = slope_ids[sig_slopes_indices]
    slopes = slopes[sig_slopes_indices]
    slope_errors = slope_errors[sig_slopes_indices]
    
    good_slope_ids = np.array([])
    
    for slope_id, slope, error in zip(slope_ids, slopes, slope_errors): 
        for i in np.linspace(slope-error, slope+error, 100):
            if i < 65:
                good_slope_ids = np.append(good_slope_ids, slope_id)
                break
    
    qm_df = pd.read_csv(qm_results)
    qm_ids = np.array(qm_df['ID'])
    qs = np.array(qm_df['rQ'])
    ms = np.array(qm_df['rM'])
    
    #q_ids = qm_ids[np.logical_and(qs<0.87, qs>0.45)]
    q_ids = qm_ids[np.where(qs<1)]
    m_ids = qm_ids[np.where(ms<-0.25)]
        
    #qpd_ids = np.intersect1d(q_ids, m_ids)
    burster_ids = np.intersect1d(q_ids, m_ids)
    
    final_ids = np.intersect1d(burster_ids, good_slope_ids)
    
    for i in final_ids: 
        lc_df = pd.read_csv(os.path.join(lc_dir, (i+'_r.csv')))
        mjds = np.array(lc_df['mjd'])
        mags = np.array(lc_df['mag'])
                
        plt.scatter(mjds, mags, color='#408ee0', s=4)
        plt.gca().invert_yaxis()
        plt.xlabel('Date (MJD)')
        plt.ylabel('r (mag)')
        plt.title(i)
        plt.savefig(os.path.join(plot_dir,(i+'.png')))
        plt.clf()
        plt.close()    
        
    #final_slopes = slopes[np.where(slopes==final_ids)]
    #final_qs = qs[np.where(qm_ids==final_ids)]
    #final_ms = ms[np.where(qm_ids==final_ids)]
    
    print(final_ids)

lc_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/light_curves/'    
qm_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/qm_plus_results.csv'
odr_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/eta/q_work_with_owen/actual_q_work/95_odr_results.csv'
#plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/aug2-7/burster_maybes'
#find_objects(lc_dir=lc_dir, qm_results=qm_results, odr_results=odr_results, plot_dir=plot_dir)        

'''
cmd_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/cmd/'
plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/slopes_demo_plots'
slopes = [3.6602377592858373, 3.6589457101373295, 3.4387383931176716]
intercepts = [9.909792577788114, 8.465106990805705, 12.530546532882909]
#slopes_demo_plot(ids=['2MASS_J20582381+4353114','FHK_163', 'FHK_405'], lc_dir=lc_dir, cmd_dir=cmd_dir, slopes=slopes, intercepts=intercepts, plot_dir=plot_dir)
'''
cmd_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/cmd'
slopes_2 = [1.9361835024995053, 0.985578160011938]
intercepts_2 = [14.939790288192524, 14.667908951835173]
plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/aug2-7/plots_to_email_to_self'
slopes_demo_plot_2(ids=['FHK_32','FHK_142'], lc_dir=lc_dir, cmd_dir=cmd_dir, slopes=slopes_2, intercepts=intercepts_2, plot_dir=plot_dir)
