def find_outlier(key,data_temp):
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Action
    key_df = pd.read_csv(key)
    ids = np.array(key_df['ID'])
    
    mean_mags = np.array([])
    mean_stds = np.array([])
    
    for obj_id in ids:
        data_df = pd.read_csv(data_temp.replace("+++",obj_id))
        mean_mag = np.mean(np.array(data_df['mag']))
        mean_magerr = np.mean(np.array(data_df['magerr']))
        
        mean_mags = np.append(mean_mags,mean_mag)
        mean_stds = np.append(mean_stds,mean_magerr)
        
        if 18.25 > mean_mag > 18 and mean_magerr>0.06:
            print(obj_id)
            print(np.max(data_df['mag'])-np.min(data_df['mag']))
        
    # Plot
    plt.rcParams['font.family']='serif'
    plt.gca().tick_params(axis='both',which='both',direction='in')
    plt.minorticks_on()
    
    plt.scatter(mean_mags,mean_stds,color='#408ee0',edgecolor='black')
    plt.gca().invert_xaxis()
    plt.ylabel(r"$\sigma_{r}$")
    plt.xlabel(r"$r$"+' (mag)')
    plt.show()
        
        
                
        
    
    
key = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/key.csv'
data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
find_outlier(key,data_temp)