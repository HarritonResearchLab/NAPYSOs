import numpy as np
from astropy.convolution import Box1DKernel, convolve
from astropy.timeseries import LombScargle
from pandas import Series, concat, date_range, read_csv, to_datetime
from scipy.interpolate import interp1d

def best_per(JD, mag, err):
    #Import(s)
    from astropy.timeseries import LombScargle
    from pandas import Series, concat, date_range, read_csv, to_datetime
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
        pass

def quas_per(JD, mag, err, ref_file):
    #Import(s)
    import pandas as pd
    from scipy.optimize import curve_fit
    import numpy as np

    

    #Action
    # Find the best period for the lc
    per = best_per(JD, mag, err)
    
    if per != 'NaN':
        
        '''SIG CALCULATION AREA'''
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
        '''SIG CALCULATION AREA'''

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

        return quas_per
    else: 
        return 'NaN'