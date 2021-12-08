def process_data_for_it(raw_lc_dir, output_dir):
    # Import(s)
    import os
    import numpy as np
    import pandas as pd
    from astropy.io import ascii
    
    # Definitions
    
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
    
    # Action
    r_tbl = os.path.join(raw_lc_dir, 'FHK_176_r.tbl')
    g_tbl = os.path.join(raw_lc_dir, 'FHK_176_g.tbl')
    raw_r_data = ascii.read(r_tbl, format='ipac', delimiter='|')
    raw_g_data = ascii.read(r_tbl, format='ipac', delimiter='|')
    
    raw_r_mjds = np.array(raw_r_data['mjd'])
    raw_r_mags = np.array(raw_r_data['mag'])
    raw_r_magerrs = np.array(raw_r_data['magerr'])
    
    raw_g_mjds = np.array(raw_g_data['mjd'])
    raw_g_mags = np.array(raw_g_data['mag'])
    raw_g_magerrs = np.array(raw_g_data['magerr'])
    
    # Apply data cuts
    
    no_streaks_mask = np.logical_or(raw_r_mjds<58448.0,raw_r_mjds>58456.0)
    raw_r_mjds = raw_r_mjds[no_streaks_mask]
    raw_r_mags = raw_r_mags[no_streaks_mask]
    raw_r_magerrs = raw_r_magerrs[no_streaks_mask]
    
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
    r_mask = np.logical_and((raw_r_mags>mean_raw_r-5*std_raw_r),(raw_r_mags<mean_raw_r+5*std_raw_r))
    r_mjds = raw_r_mjds[r_mask]
    r_mags = raw_r_mags[r_mask]
    r_magerrs = raw_r_magerrs[r_mask]

    g_mask = np.logical_and((raw_g_mags>mean_raw_g-5*std_raw_g),(raw_g_mags<mean_raw_g+5*std_raw_g))
    g_mjds = raw_g_mjds[g_mask]
    g_mags = raw_g_mags[g_mask]
    g_magerrs = raw_g_magerrs[g_mask]
    
    print('NUM OBS: ', len(r_mags))
    print('MEAN R: ', np.mean(r_mags))    
    print('SOKOLOVSKY NU: ', sokolovskyNu(r_mags, r_magerrs))
    
    r_zipped = list(zip(r_mjds, r_mags, r_magerrs))
    r_df = pd.DataFrame(r_zipped, columns=['mjd', 'mag', 'magerr'])
    r_df.to_csv(output_dir+'/FHK_176_r.csv', index=False)
    g_zipped = list(zip(g_mjds, g_mags, g_magerrs))
    g_df = pd.DataFrame(g_zipped, columns=['mjd', 'mag', 'magerr'])
    g_df.to_csv(output_dir+'/FHK_176_g.csv', index=False)
    
raw_lc_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves' 
output_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/light_curves'
process_data_for_it(raw_lc_dir, output_dir)
