from numpy.core.numeric import full


def make_plots(full_dir, final_dir):
    # Import(s)
    import os
    import pandas as pd 
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    from scipy import interpolate
    
    # Action
    
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
    
    full_samp_mean_mags = np.array([])
    full_samp_mean_magerrs = np.array([])
    full_samp_nus = np.array([])
    variable_mean_mags = np.array([])
    variable_mean_magerrs = np.array([])
    variable_nus = np.array([])

    for full_samp_file in os.listdir(full_dir):
        if '_r.csv' in full_samp_file:
            data_file = pd.read_csv(os.path.join(full_dir, full_samp_file))
            mags = data_file['mag']
            magerrs = data_file['magerr']
            full_samp_mean_mags = np.append(full_samp_mean_mags, 
                                            np.mean(mags))
            full_samp_mean_magerrs = np.append(full_samp_mean_magerrs, 
                                               np.mean(magerrs))
            full_samp_nus = np.append(full_samp_nus, sokolovskyNu(mags,
                                                                  magerrs))
    
    for full_samp_file in os.listdir(final_dir):
        if '_r.csv' in full_samp_file:
            data_file = pd.read_csv(os.path.join(full_dir, full_samp_file))
            mags = data_file['mag']
            magerrs = data_file['magerr']
            variable_mean_mags = np.append(variable_mean_mags, 
                                            np.mean(mags))
            variable_mean_magerrs = np.append(variable_mean_magerrs, 
                                               np.mean(magerrs))
            variable_nus = np.append(variable_nus, sokolovskyNu(mags,
                                                                  magerrs))                                        
    
    
    # mean mags vs mean magerrs
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    plt.rcParams['font.family'] = 'serif'
    
    plt.scatter(full_samp_mean_mags,full_samp_mean_magerrs,color='lightsteelblue',edgecolors='black', label='Initial Sample')
    
    plt.scatter(variable_mean_mags,variable_mean_magerrs,color='#408ee0',edgecolors='black', label='Final Sample')
    plt.tick_params(axis='both',which='both',direction='in')
    plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
    plt.legend(loc='upper right',fontsize=7,fancybox=False,edgecolor='black',shadow=False)
    plt.gca().invert_xaxis()
    plt.xlabel(r'$r$'+' (mag)')
    plt.ylabel(r'$\sigma_{r}$')
    plt.savefig('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/plots for lah edits/mag_magerr.png',dpi=500)
    plt.clf()
    plt.close()
    
    # mean mag vs nu 
    inter_df = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/inter_bounds.csv')
    inter_x = np.array(inter_df['mean mag'])
    inter_y = np.array(inter_df['nu cutoff'])
    
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    plt.xlabel(r'$r$'+' (mag)')
    plt.ylabel('Peak-to-peak variability ('+r'$\nu$'+')')
    plt.scatter(full_samp_mean_mags,full_samp_nus,color='lightsteelblue',edgecolors='black', label='Initial Sample')
    plt.scatter(variable_mean_mags,variable_nus,color='#408ee0',edgecolors='black', label='Final Sample')
    plt.plot(inter_x,inter_y,color='black',lw=1,ls='--',label='15th Percentile')
    plt.yscale('log')
    plt.xlim(0.99*min(variable_mean_mags),1.01*max(variable_mean_mags))
    plt.gca().invert_xaxis()
    plt.legend(loc='lower left',fontsize=7,fancybox=False,edgecolor='black',shadow=False)
    plt.tick_params(axis='both',which='both',direction='in')
    plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
    plt.savefig('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/plots for lah edits/nu_cutoff.png',dpi=500)
    plt.clf()
    plt.close()
    
full_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves'
final_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves'
make_plots(full_dir, final_dir)