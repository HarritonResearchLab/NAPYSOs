import numpy as np 
import pandas as pd
import uncertainties as u  
from uncertainties import unumpy
from progress.bar import Bar



key_df = pd.read_csv(r'C:\Users\Research\Documents\GitHub\NAPYSOs\recovered\LAH2.0\data\AUGUST_5th\key.csv')
ids = np.array(key_df['ID'])

bar = Bar('Processing', max=len(ids))

def interpolate(x, x1, y1, x2, y2): 
    m = (y2-y1)/(x2-x1)
    b = y1-m*x1
    return x*m+b

for id in ids: 
    r_path = r'C:\Users\Research\Documents\GitHub\NAPYSOs\recovered\LAH2.0\data\AUGUST_5th\light_curves\***_r.csv'.replace('***',id)
    r_df = pd.read_csv(r_path)

    g_path = r'C:\Users\Research\Documents\GitHub\NAPYSOs\recovered\LAH2.0\data\AUGUST_5th\light_curves\***_g.csv'.replace('***',id)
    g_df = pd.read_csv(g_path)

    g_mjds = np.array(g_df['mjd'])
    g_mags = np.array(g_df['mag'])
    g_magerrs = np.array(g_df['magerr'])

    r_mjds = np.array(r_df['mjd'])
    r_mags = np.array(r_df['mag'])
    r_magerrs = np.array(r_df['magerr'])

    r_idx_to_remove = np.logical_or(r_mjds>np.max(g_mjds), r_mjds<np.min(g_mjds))
    r_mjds, r_mags, r_magerrs = (np.delete(arr,r_idx_to_remove) for arr in [r_mjds, r_mags, r_magerrs])

    interpolated_mags = np.array([])
    interpolated_magerrs = np.array([])
    interpolated_dates = np.array([])
    colors = np.array([])
    color_errs = np.array([])

    for r_index, r_mjd in enumerate(r_mjds): 
        diffs = g_mjds - r_mjd
        lowers = diffs[np.where(diffs<0)] 
        uppers = diffs[np.where(diffs>0)] 
        
        if len(lowers) == 0: 
            lower = 0
        
        else: 
            lower = np.max(lowers)
        
        if len(uppers) == 0: 
            upper = 0
        
        else: 
            upper = np.min(uppers)

        window = np.abs(upper+lower)  

        if window < 3:
            if lower+upper == 0: 
                print(lower, upper)
            
            lower_index = np.where(diffs==lower)
            upper_index = np.where(diffs==upper)

            r_mjd_u = u.ufloat(r_mjd, 0)
            g_mjd_1 = u.ufloat(g_mjds[lower_index], 0)
            g_mag_1 = u.ufloat(g_mags[lower_index], g_magerrs[lower_index])
            g_mjd_2 = u.ufloat(g_mjds[upper_index], 0)
            g_mag_2 = u.ufloat(g_mags[upper_index], g_magerrs[upper_index])
            
            interpolated_mag_u = interpolate(r_mjd_u, g_mjd_1, g_mag_1, g_mjd_2, g_mag_2)
            interpolated_mag = u.nominal_value(interpolated_mag_u)
            interpolated_magerr = u.std_dev(interpolated_mag_u)
            
            color = interpolated_mag_u - u.ufloat(r_mags[r_index], r_magerrs[r_index])

            interpolated_mags = np.append(interpolated_mags, interpolated_mag)
            interpolated_magerrs = np.append(interpolated_magerrs, interpolated_magerr)
            interpolated_dates = np.append(interpolated_dates, r_mjd)
            colors = np.append(colors, u.nominal_value(color))
            color_errs = np.append(color_errs, u.std_dev(color))

    zipped = list(zip(interpolated_dates, interpolated_mags, interpolated_magerrs, 
                      colors, color_errs, r_mags, r_magerrs))

    cmd_df_columns = ['date (mjd)','g_mag','g_magerr','color','color_err','r_mag','r_magerr']

    cmd_df = pd.DataFrame(zipped, columns=cmd_df_columns)
    
    cmd_path = r"C:\Users\Research\Documents\GitHub\NAPYSOs\referee_changes\color-mag\reprocessing\cmd_files"
    cmd_path += r"\\" + id + "_cmd.csv" 
    cmd_df.to_csv(cmd_path, index=False)
    bar.next()

bar.finish()




