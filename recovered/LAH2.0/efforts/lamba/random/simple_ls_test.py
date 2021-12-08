from re import A


def simple_ls(id,data_temp): 
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from astropy.timeseries import LombScargle
    
    # Action
    data_df = pd.read_csv(data_temp.replace('+++',id))
    mjds = np.array(data_df['mjd'])
    mags = np.array(data_df['mag'])
    
    minf = 1/250
    maxf = 1/0.5

    # Periodogram of light curve
    periodogram = LombScargle(mjds,mags)
    ls_freqs,ls_powers = periodogram.autopower(method='fast',minimum_frequency=minf,maximum_frequency=maxf)
    ls_periods = 1/ls_freqs
    
    best_index = np.argmax(ls_powers)
    best_per = float(ls_periods[best_index])
    best_pow = float(ls_powers[best_index])
    
    print(best_per)
    
    phased_dates = (np.mod(mjds, best_per))/best_per 
    phased_dates_cycle_2 = phased_dates + 1
    
    # Plot
    plt.style.use('seaborn-darkgrid')
    plt.rcParams['font.family']='serif'
    fig, axs = plt.subplots(2,2)
    
    axs[0,0].plot(ls_periods,ls_powers,label='Best Period: '+str(round(best_per,2))+'\nBest Bhardwaj: 3.95',color='#408ee0')
    axs[0,0].legend(fontsize='7')
    axs[0,0].set_xlabel('Period (d)')
    axs[0,0].set_ylabel('Normalized Power')
    axs[0,0].set_xscale('log')
    
    axs[0,1].scatter(phased_dates,mags,color='#408ee0',s=5)
    axs[0,1].scatter(phased_dates_cycle_2,mags,color='#408ee0',s=5)
    axs[0,1].set_ylabel('Mag')
    axs[0,1].set_xlabel('Cycle')
    axs[0,1].invert_yaxis()

    axs[1,0].axis('off')
    axs[1,1].axis('off')
    
    plt.subplots_adjust(wspace=0.5)
    
    plt.show()
    
    
id = 'FHK_8'
data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
simple_ls(id,data_temp)