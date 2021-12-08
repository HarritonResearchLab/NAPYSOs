def simple_lomb_scarg(data_file,id,min_per,max_per):
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from astropy.timeseries import LombScargle
    from matplotlib.ticker import AutoMinorLocator

    # Action

    df = pd.read_csv(data_file)
    
    mjd = np.array(df['mjd'])
    mag = np.array(df['mag'])
    magerr = np.array(df['magerr'])

    # Initializing these here to optimize speed
    harmonics = np.array(5*[1])/range(1,6)
    lower_harmonics = 0.98*harmonics
    upper_harmonics = 1.02*harmonics
    lph = 1/upper_harmonics # lph = Lower Period Harmonics (for lower bounds of harmonic period ranges)
    uph = 1/lower_harmonics # uph = Upper Period Harmonics (for upper bounds of harmonic period ranges)

    def per_is_good(period):
        harmonic = False
        for harmonic_idx in range(0,len(harmonics)):
            if lph[harmonic_idx] < period < uph[harmonic_idx]:
                harmonic = True
                break
        if harmonic==True: 
            return False
        else: 
            return True

    minf = 1.0/max_per 
    maxf = 1.0/min_per
    periodogram = LombScargle(mjd,mag,magerr)
    ls_freqs,ls_powers = periodogram.autopower(method='fast',minimum_frequency=minf,maximum_frequency=maxf)
    faps = periodogram.false_alarm_level([0.1,0.05,0.01])
    
    # CLEANING
    for f in ls_freqs: 
        per = 1.0/f
        half_per = 0.5*per
        double_per = 2.0*per
        f_index = np.where(ls_freqs==f)

        if per_is_good(per)==True:
            if per_is_good(half_per)==False or per_is_good(double_per)==False:
                ls_freqs = np.delete(ls_freqs,f_index)
                ls_powers = np.delete(ls_powers,f_index)
        else: 
            ls_freqs = np.delete(ls_freqs,f_index)
            ls_powers = np.delete(ls_powers,f_index)
    
    best_index = np.argmax(ls_powers)
    best_per = 1/ls_freqs[best_index]
    best_power = ls_powers[best_index]
    periods = 1/ls_freqs


    # Fold the light curve
    phased_dates = (np.mod(mjd, best_per))/best_per 
    phased_dates_cycle_2 = phased_dates + 1

    # Plot
    plt.rcParams['font.family']='serif'
    fig, axs = plt.subplots(2,1)
    plt.suptitle(id,fontsize='medium')

    # Folded light curve plot
    xlabel = 'Phase (P = '+str(round(best_per, 3))+' d)'
    axs[0].set_xlabel(xlabel)
    axs[0].set_ylabel('Mag')
    axs[0].scatter(phased_dates,mag,color='#408ee0',s=2)
    axs[0].scatter(phased_dates_cycle_2,mag,color='#408ee0',s=2)

    axs[0].invert_yaxis()
    axs[0].set_title("Folded Light Curve", fontsize=10)
    axs[0].locator_params(axis='y', nbins=5)
    axs[0].tick_params(axis='both',which='both',direction='in')
    axs[0].tick_params(which='both',bottom=True,top=True,left=True,right=True)
    axs[0].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    axs[0].xaxis.set_minor_locator(AutoMinorLocator())
    axs[0].yaxis.set_minor_locator(AutoMinorLocator())

    # Powers plot
    axs[1].plot(periods,ls_powers,color='#408ee0',label='Power')
    for i in range(0,len(harmonics)):
        axs[1].axvspan(lph[i],uph[i],color='red',alpha=0.2)
    axs[1].set_xscale('log')
    axs[1].set_xlabel('Period d')
    axs[1].set_xscale('log')
    axs[1].set_ylabel('Power')
    axs[1].set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
    axs[1].set_ylim(0,1)
    axs[1].set_xlim(0.5,250)
    axs[1].set_title('Periodogram', fontsize=10)
    axs[1].legend(loc='upper right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
    axs[1].tick_params(axis='both',which='both',direction='in')
    axs[1].tick_params(which='both',bottom=True,top=True,left=True,right=True)
    axs[1].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    axs[1].yaxis.set_minor_locator(AutoMinorLocator())
    for i in faps: 
        axs[1].axhline(y=i,color='brown')
    


    plt.subplots_adjust(hspace=0.4)
    
    plt.show()


data_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/FHK_388_r.csv'
simple_lomb_scarg(data_file=data_file,id='FHK_388',min_per=0.5,max_per=250)