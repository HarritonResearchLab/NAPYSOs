def generate_final_sample(in_key,raw_g_temp,raw_r_temp,out_key,sigma,nu_percentile,max_mean_g,max_mean_r,min_obs,out_data_temp):
    # Import(s)
    import numpy as np
    import pandas as pd
    
    # Action 

    starting_key = pd.read_csv(in_key) 

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

    def first_three_cuts(starting_key):
        # Import(s)
        from astropy.io import ascii

        # Action
        starting_ids = list(starting_key['ID'])
        initial_sample_size = len(starting_ids)
        ids_that_meet_mag_cut = []
        more_than_thirty_obs_ids = []

        for id in starting_ids: 
            raw_g_data = ascii.read(raw_g_temp.replace('+++',id),format='ipac',delimiter='|')
            raw_g_mjds = np.array(raw_g_data['mjd'])
            raw_g_mags = np.array(raw_g_data['mag'])
            raw_g_magerrs = np.array(raw_g_data['magerr'])

            raw_r_data = ascii.read(raw_r_temp.replace('+++',id),format='ipac',delimiter='|')
            raw_r_mjds = np.array(raw_r_data['mjd'])
            raw_r_mags = np.array(raw_r_data['mag'])
            raw_r_magerrs = np.array(raw_r_data['magerr'])

            if len(raw_g_mjds)>5 and len(raw_r_mjds)>5: 
                mean_raw_g = np.mean(raw_g_mags)
                std_raw_g = np.std(raw_g_mags)

                mean_raw_r = np.mean(raw_r_mags)
                std_raw_r = np.std(raw_r_mags)

                # Ignore streaks range from r files
                no_streaks_mask = np.logical_or(raw_r_mjds<58448.0,raw_r_mjds>58456.0)
                raw_r_mjds = raw_r_mjds[no_streaks_mask]
                raw_r_mags = raw_r_mags[no_streaks_mask]
                raw_r_magerrs = raw_r_magerrs[no_streaks_mask]

                # Clip r and g files
                r_mask = np.logical_and((raw_r_mags>mean_raw_r-sigma*std_raw_r),(raw_r_mags<mean_raw_r+sigma*std_raw_r))
                r_mjds = raw_r_mjds[r_mask]
                r_mags = raw_r_mags[r_mask]
                r_magerrs = raw_r_magerrs[r_mask]

                g_mask = np.logical_and((raw_g_mags>mean_raw_g-sigma*std_raw_g),(raw_g_mags<mean_raw_g+sigma*std_raw_g))
                g_mjds = raw_g_mjds[g_mask]
                g_mags = raw_g_mags[g_mask]
                g_magerrs = raw_g_magerrs[g_mask]

                if len(r_mags)>5 and len(g_mags) > 5:
                    # Test if good
                    mean_g = np.mean(g_mags)
                    mean_r = np.mean(r_mags)
                    
                    if mean_g < max_mean_g and mean_r < max_mean_r:
                        ids_that_meet_mag_cut.append(id)

                        if len(g_mjds) > min_obs and len(r_mjds) > min_obs: 
                            more_than_thirty_obs_ids.append(id)
                            
                            # Save cleaned data file 
                            g_zipped = list(zip(g_mjds,g_mags,g_magerrs))
                            r_zipped = list(zip(r_mjds,r_mags,r_magerrs))
                            col_names =['mjd','mag','magerr']

                            g_temp_out_df = pd.DataFrame(data=g_zipped,columns=col_names)
                            g_temp_out_df.to_csv(out_data_temp.replace('+++',(id+'_g')),index=False)
                            r_temp_out_df = pd.DataFrame(data=r_zipped,columns=col_names)
                            r_temp_out_df.to_csv(out_data_temp.replace('+++',(id+'_r')),index=False)

        # Print results summary
        print('\n'+'Results summary:'+'\n')
        print('Criteria      |      Remaining Sample size')
        print('----------------------------------------')
        print('Initial       |      '+str(initial_sample_size))
        print('----------------------------------------')
        print('Mag sigma cut |      '+str(len(ids_that_meet_mag_cut)))
        print('----------------------------------------')
        print('Num obs cut   |      '+str(len(more_than_thirty_obs_ids)))

        # Return ids that made first three cuts
        return more_than_thirty_obs_ids

    almost_ready_ids = first_three_cuts(starting_key=starting_key)

    def nu_variability_cut(nu_percentile):
        # Import(s)
        from scipy import interpolate
        import os
        import matplotlib.pyplot as plt
        from matplotlib.ticker import AutoMinorLocator

        # Action 
        final_ids = []
        mean_mags = []
        mean_magerrs = []
        nus = []
        

        # Set mag vs nu plot up
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['mathtext.fontset'] = 'dejavuserif'
        plt.xlabel(r'$r$'+' (mag)')
        plt.ylabel('Peak-to-peak variability ('+r'$\nu$'+')')
        # Continue

        for id in almost_ready_ids: 
            data_file = pd.read_csv(out_data_temp.replace('+++',(id+'_r')))
            mags = np.array(data_file['mag'])
            magerrs = np.array(data_file['magerr'])
            mean_mag = np.mean(mags)
            nu = sokolovskyNu(x=mags,xerr=magerrs)
            mean_mags.append(mean_mag)
            nus.append(nu)
            mean_magerrs.append(float(np.mean(magerrs)))

        min_mean = int(np.floor(min(mean_mags)))
        max_mean = int(np.ceil(max(mean_mags)))

        lower_bounds = list(range(min_mean,max_mean))
        median_bounds = []
        cutoffs = []

        mean_mags = np.array(mean_mags)
        nus = np.array(nus)

        # Continue plot
        plt.scatter(mean_mags,nus,color='#408ee0',edgecolors='black',label='NAP YSOs')

        # Continue with the rest

        for lower_bound in lower_bounds: 
            upper_bound = lower_bound + 1.0
            mask = np.logical_and(mean_mags<upper_bound,mean_mags>=lower_bound)
            masked_nus = nus[mask]
            cutoff = float(np.percentile(masked_nus,q=[nu_percentile]))
            cutoffs.append(cutoff)
            
            med_bound = np.median([lower_bound,upper_bound])
            median_bounds.append(med_bound)    
        
        
        f = interpolate.interp1d(median_bounds,cutoffs,kind='linear',fill_value='extrapolate')
        inter_x = np.linspace(min(mean_mags),max(mean_mags),50)
        inter_y = f(inter_x)
        
    
        ziPPeD = list(zip(inter_x, inter_y))    
        for_aug_df = pd.DataFrame(ziPPeD, columns=['mean mag', 'nu cutoff'])
        for_aug_df.to_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/inter_bounds.csv', index=False)
        
        # Continue with plot
        plt.plot(inter_x,inter_y,color='black',lw=1,ls='--',label='15th Percentile')
        plt.yscale('log')
        plt.xlim(0.99*min(mean_mags),1.01*max(mean_mags))
        plt.gca().invert_xaxis()
        plt.legend(loc='lower left',fontsize=7,fancybox=False,edgecolor='black',shadow=False)
        plt.tick_params(axis='both',which='both',direction='in')
        plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
        plt.savefig('/home/thaddaeus/nu_cutoff.png',dpi=500)
        plt.clf()
        plt.close()

        

        # Continue with the rest
        variable_mean_mags = []
        variable_mean_magerrs = []
        for id, mean_mag, nu, magerr in zip(almost_ready_ids,list(mean_mags),list(nus),mean_magerrs):
            if nu > f(mean_mag):
                final_ids.append(id)
                variable_mean_mags.append(mean_mag)
                variable_mean_magerrs.append(magerr)
            else: 
                os.remove(out_data_temp.replace('+++',(id+'_r')))
                os.remove(out_data_temp.replace('+++',(id+'_g')))

        # Mean mags vs mean magerrs plot
        plt.rcParams['mathtext.fontset'] = 'dejavuserif'
        plt.rcParams['font.family'] = 'serif'
        plt.scatter(variable_mean_mags,variable_mean_magerrs,color='#408ee0',edgecolors='black')
        plt.tick_params(axis='both',which='both',direction='in')
        plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
        plt.gca().invert_xaxis()
        plt.xlabel(r'$r$'+' (mag)')
        plt.ylabel(r'$\sigma_{r}$')
        plt.savefig('/home/thaddaeus/mag_magerr.png',dpi=500)

        # Report results
        print('----------------------------------------')
        print('Nu            |      '+str(len(final_ids)))

        # Return final id list
        return final_ids

    final_ids = nu_variability_cut(nu_percentile=nu_percentile)
    
    # Get RAs and DECs for final ids
    final_ras = []
    final_decs = []

    all_ids = np.array(starting_key['ID'])
    all_ras = np.array(starting_key['RA'])
    all_decs = np.array(starting_key['DEC'])
    
    for id in final_ids: 
        all_ids_index = np.where(all_ids==id)
        ra = all_ras[all_ids_index]
        ra = float(ra[0])
        dec = all_decs[all_ids_index]
        dec = float(dec[0])
        final_ras.append(str(ra))
        final_decs.append(str(dec))

    final_df = pd.DataFrame(data=list(zip(final_ids,final_ras,final_decs)),columns=['ID','RA','DEC'])
    final_df.to_csv(out_key,index=False)




i_key='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/key.csv'
raw_g_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/raw_from_lah/ZTF_G/+++_g.tbl'
raw_r_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/raw_from_lah/ZTF_R/+++_r.tbl'
out_data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/light_curves/+++.csv'
out_key = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/key.csv'


import time
start = time.time()
generate_final_sample(in_key=i_key,raw_g_temp=raw_g_temp,raw_r_temp=raw_r_temp,out_key=out_key,sigma=5,nu_percentile=15.0,max_mean_g=20.8,max_mean_r=20.6,min_obs=30,out_data_temp=out_data_temp)
end = time.time()
print('\nTime to process: '+str(round((end-start)/60,4))+' minutes.')