def print_bounds(key,r_data_temp,g_data_temp): 
    # Import(s)
    from astropy.io import ascii
    import numpy as np
    import pandas as pd
    
    # Action
    key_df = pd.read_csv(key)
    ids = np.array(key_df["ID"])
    
    mins = np.array([])
    maxs = np.array([])
    for id in ids: 
        r_table = ascii.read(r_data_temp.replace('+++',id),format='ipac',delimiter='|')
        r_dates = np.array(r_table['mjd'])
        
        if len(r_dates>2):
            r_min = np.amin(r_dates)
            r_max = np.amax(r_dates)
            mins = np.append(mins,r_min)
            maxs = np.append(maxs,r_max)
        
        g_table = ascii.read(g_data_temp.replace('+++',id),format='ipac',delimiter='|')
        g_dates = np.array(g_table['mjd'])
        if len(g_dates)>2:
            g_min = np.amin(g_dates)
            g_max = np.amax(g_dates)
            mins = np.append(mins,g_min)
            maxs = np.append(maxs,g_max)
        
    print(np.amin(mins),np.amax(maxs))
            

key = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/key.csv'
r_data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/raw_from_lah/ZTF_R/+++_r.tbl'
g_data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/raw_from_lah/ZTF_G/+++_g.tbl'
print_bounds(key,r_data_temp,g_data_temp)