import numpy as np
from astropy.convolution import Box1DKernel, convolve
from astropy.timeseries import LombScargle
from pandas import Series, concat, date_range, read_csv, to_datetime
from scipy.interpolate import interp1d

def codyM(x):
    #Import(s)
    from scipy import stats 
    import numpy as np
    
    #Action
    x = np.array(x)
    m_metric = (np.mean([stats.mstats.mquantiles(x,prob=0.9),stats.mstats.mquantiles(x,prob=0.1)])-np.median(x))/np.sqrt(((x-x.mean())**2).sum()/len(x))
    return m_metric

def best_per(JD, mag, err):
    #Import(s)
    from astropy.timeseries import LombScargle
    from pandas import Series, concat, date_range, to_datetime
    import numpy as np

    #Action

    # Generate a Dirac Comb, our window function
    time = to_datetime(JD, unit="D", origin="julian")
    time_not_obs = date_range(time.min(), time.max(), periods=1000)
    base = Series(np.zeros(len(time_not_obs)), index=time_not_obs)
    teeth = Series(np.ones(len(time)), index=time)
    dirac_comb = concat([base, teeth]).sort_index()

    # We only want to look for signals that we have 3 periods observed. ###I SHOULD ADD THIS TO OUR LS FUNC
    maxp = (JD.max() - JD.min()) / 3
    minf = 1/maxp

    # First, we find the periodigram of the window function
    JD_W = dirac_comb.index.to_julian_date()
    mag_W = dirac_comb.values
    periodogram_W = LombScargle(JD_W, mag_W)
    freq_W, power_W = periodogram_W.autopower(minimum_frequency=minf, maximum_frequency=1/0.1)

    # Now, for the original lightcurve
    periodogram = LombScargle(JD, mag)
    freq, power = periodogram.autopower(minimum_frequency=minf, maximum_frequency=1/0.1)

    # Mask out peak window-function frequencies from the data with a notch
    # width of 0.03 Hz on either side.
    high_power_W = power_W.mean() + 2 * power_W.std()
    for idx in np.argwhere(power_W > high_power_W):
        v = freq[idx]
        dv = 0.03
        vmin, vmax = v - dv, v + dv
        power[np.logical_and(freq > vmin, freq < vmax)] = 0

    # We find the best frequency that has a false alarm probability < 0.1.
    faps = periodogram.false_alarm_probability(power)
    try: 
        idx = np.argmax(power)
        best_fap = faps[idx]
        while best_fap > 0.1:
            power = np.delete(power, idx)
            idx = np.argmax(power)
            best_fap = faps[idx]
        best_freq = freq[idx]

        # Periods are more human-friendly :)
        return 1 / best_freq
    except ValueError:
        return 'NaN'
        

def quas_per(JD, mag, err, ref_file):
    #Import(s)
    import pandas as pd
    from scipy.optimize import curve_fit
    import numpy as np
    from astropy.convolution import Box1DKernel, convolve
    

    #Action
    # Find the best period for the lc
    per = best_per(JD, mag, err)
    
    if per != 'NaN':

        #Janky way to get sig(?)
        ref_df = pd.read_csv(ref_file)
        mean_mags = np.array(ref_df['mean_mag'])
        mean_stds = np.array(ref_df['mean_std'])
        
        def exponential(x,a,b,c):
            return((a*np.exp(b*x))+c)
        
        pars, cov = curve_fit(f=exponential,xdata=mean_mags,ydata=mean_stds)
        
        def bf_exp(x): 
            return exponential(x,pars[0],pars[1],pars[2])
        
        sig = bf_exp(np.mean(mag))
        

        # Shifting gears entirely, we now must create the residual curve
        phase = JD % per
        mag = mag[np.argsort(phase)]

        # We use three periods and extract the middle to prevent edge effects
        three_periods = np.concatenate((mag, mag, mag))
        boxcar = Box1DKernel(len(mag) // 4)
        smooth_mag = convolve(three_periods, boxcar)
        smooth_mag = smooth_mag[np.size(mag):2*np.size(mag)]

        resid_mag = mag - smooth_mag

        quas_per = (
            (np.nanstd(resid_mag) ** 2 - sig ** 2)
            / (np.nanstd(mag) - sig ** 2)
        )

        #return quas_per, resid_mag, smooth_mag, per
        return quas_per
    else: 
        return 'NaN'

def create_sig(JD, mag, err, ref_file):
    #Import(s)
    import pandas as pd
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import numpy as np

    #Action

    df = pd.read_csv(ref_file)
    mean_mags = np.array(df['mean_mag'])
    mean_stds = np.array(df['mean_std'])

    def exponential(x,a,b,c):
        return((a*np.exp(b*x))+c)
    
    pars, cov = curve_fit(f=exponential,xdata=mean_mags,ydata=mean_stds)
    
    def bf_exp(x): 
        return exponential(x,pars[0],pars[1],pars[2])

    estimated_y = bf_exp(mean_mags)

    plt.scatter(mean_mags,mean_stds,color='cornflowerblue',label='Real data')
    plt.scatter(mean_mags,estimated_y,color='orange',label='Exponential best fit')

    plt.xlabel('Mean mag')
    plt.ylabel('Mean std')
    plt.legend()
    plt.show()
    
    sig = bf_exp(np.mean(mag))
    
    return sig

def create_ref_file():
    #Import(s)
    import numpy as np
    import pandas as pd

    #Action
    mean_mags = []
    mean_stds = []
    ids = list(pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/key.csv')['ID'])
    for id in ids:
        data_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++_r.csv'.replace('+++',id)
        df = pd.read_csv(data_file)
        r_mags = np.array(df['mag'])
        r_magerrs = np.array(df['magerr'])

        mean_mags.append(float(np.mean(r_mags)))
        mean_stds.append(float(np.mean(r_magerrs)))

    #Sort them

    sorted_mean_mags = sorted(mean_mags)
    sorted_mean_stds = []
    for mean_mag in sorted_mean_mags:
        item_index = mean_mags.index(mean_mag)
        sorted_mean_stds.append(mean_stds[item_index])
    
    zipped_list = list(zip(sorted_mean_mags,sorted_mean_stds))
    out_df = pd.DataFrame(data=zipped_list,columns=['mean_mag','mean_std'])
    out_df.to_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/mean_mags_and_stds.csv',index=False)
        
def figure_2(key,ref_file):
    #Import(s)
    import pandas as pd
    import matplotlib.pyplot as plt
    from progress.bar import Bar
    from shutil import copyfile

    ms = []
    qs = []

    #Action
    ids = list(pd.read_csv(key)['ID'])
    bar = Bar('Processing...',max=len(ids))
    
    ids_with_metrics = []
    for id in ids: 
        r_file = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++_r.csv'.replace('+++',id))
        r_mjds = np.array(r_file['mjd'])
        JD = np.array(r_mjds)+2400000.5
        r_mags = np.array(r_file['mag'])
        r_magerrs = np.array(r_file['magerr'])
        q = quas_per(JD=JD,mag=r_mags,err=r_magerrs,ref_file=ref_file)
        if q != 'NaN':
            m = codyM(r_mags)
            ms.append(m)
            qs.append(q)
            ids_with_metrics.append(id)

        bar.next()

    bar.finish()
    
    zipped = list(zip(ids_with_metrics,qs,ms)) 
    m_vs_q_df = pd.DataFrame(data=zipped,columns=['ID','Q','M'])
    m_vs_q_df.to_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/gamma/testing_qm/q_vs_m_results.csv',index=False)

    #Plot
    plt.rcParams['font.family'] = 'serif'
    plt.scatter(qs,ms,marker='.',color='cornflowerblue',s=2)
    plt.axhline(y=0.25,xmin=0,xmax=1,linewidth=1,color='indianred',linestyle='--')
    plt.axhline(y=-0.25,xmin=0,xmax=1,linewidth=1,color='indianred',linestyle='--')
    plt.axvline(x=0.15,ymin=-1,ymax=1,linewidth=1,color='indianred',linestyle='--')
    plt.axvline(x=0.85,ymin=-1,ymax=1,linewidth=1,color='indianred',linestyle='--')
    plt.xlim(0,1)
    plt.ylim(-1,1)
    plt.gca().invert_yaxis()
    plt.xlabel('Quasi-Periodicity (Q)')
    plt.ylabel('Flux Asymmetry (M)')
    plt.show()

def figure_2_plot_only(ref_file):
    #Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from progress.bar import Bar

    #Action

    df = pd.read_csv(ref_file)
    qs = list(np.array(df['Q'])*1.3)
    ms = list(df['M'])

    #Plot
    plt.rcParams['font.family'] = 'serif'
    plt.scatter(qs,ms,marker='.',color='cornflowerblue',s=2)
    plt.axhline(y=0.25,xmin=0,xmax=1,linewidth=1,color='indianred',linestyle='--')
    plt.axhline(y=-0.25,xmin=0,xmax=1,linewidth=1,color='indianred',linestyle='--')
    plt.axvline(x=0.15,ymin=-1,ymax=1,linewidth=1,color='indianred',linestyle='--')
    plt.axvline(x=0.85,ymin=-1,ymax=1,linewidth=1,color='indianred',linestyle='--')
    plt.xlim(0,1)
    plt.ylim(-1,1)
    plt.gca().invert_yaxis()
    plt.xlabel('Quasi-Periodicity (Q)')
    plt.ylabel('Flux Asymmetry (M)')
    plt.show()

def investigate_q(data_file):
    #Import(s)
    import pandas as pd 
    import matplotlib.pyplot as plt

    #Action

    r_file = pd.read_csv(data_file)
    r_mjds = np.array(r_file['mjd'])
    r_mags = np.array(r_file['mag'])
    r_magerrs = np.array(r_file['magerr'])
    JD = np.array(r_mjds)+2400000.5
    q = quas_per(JD=JD,mag=r_mags,err=r_magerrs,ref_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/mean_mags_and_stds.csv')

    #Plot

    fig, axs = plt.subplots(2,2,sharex=True)

    axs[0,0].scatter(r_mjds,r_mags,s=2)
    axs[0,0].set_xlabel('Date (MJD)')
    axs[0,0].set_ylabel('Mag')
    axs[0,0].invert_yaxis()

    axs[0,1].scatter(r_mjds,q[1],s=2)
    axs[0,1].set_xlabel('Date (MJD)')
    axs[0,1].set_ylabel('Residual Mag')
    axs[0,1].invert_yaxis()

    axs[1,0].scatter(r_mjds,q[2],s=2)
    axs[1,0].set_xlabel('Date (MJD)')
    axs[1,0].set_ylabel('Smoothed Mag')
    axs[1,0].invert_yaxis()

    plt.subplots_adjust(hspace=0.5,wspace=0.5)

    plt.show()

investigate_q('/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/FHK_282_r.csv')

