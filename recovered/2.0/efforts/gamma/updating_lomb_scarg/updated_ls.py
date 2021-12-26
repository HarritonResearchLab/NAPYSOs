def updated_ls(mjds,mags,magerrs,min_period,fap_levels,out_type):
    #Import(s)
    import numpy as np
    import pandas as pd
    from pandas import Series, concat, date_range, to_datetime
    from astropy.timeseries import LombScargle
    
    #Action
    
    JD = np.array(mjds)+2400000.5


    # Generate a Dirac Comb, our window function
    time = to_datetime(JD, unit="D", origin="julian")
    time_not_obs = date_range(time.min(), time.max(), 
    iods=1000)
    base = Series(np.zeros(len(time_not_obs)), index=time_not_obs)
    teeth = Series(np.ones(len(time)), index=time)
    dirac_comb = concat([base, teeth]).sort_index()

    # Max period must have at least 3 periods observed.
    maxp = (JD.max() - JD.min()) / 3
    minf = 1/maxp

    # Periodigram of the window function
    JD_W = dirac_comb.index.to_julian_date()
    mags_W = dirac_comb.values
    periodogram_W = LombScargle(JD_W, mags_W)
    freq_W, power_W = periodogram_W.autopower(method='fastchi2',minimum_frequency=minf, maximum_frequency=1/min_period)

    # Periodogram for the original lightcurve
    periodogram = LombScargle(JD,mags,magerrs)
    freq, power = periodogram.autopower(method='fastchi2',minimum_frequency=minf, maximum_frequency=1/min_period)

    # Mask out peak window-function frequencies from the data with a notch
    # width of 0.03 Hz on either side.
    high_power_W = power_W.mean() + 2 * power_W.std()
    for idx in np.argwhere(power_W > high_power_W):
        v = freq[idx]
        dv = 0.03
        vmin, vmax = v - dv, v + dv
        power[np.logical_and(freq > vmin, freq < vmax)] = 0

    # We find the best frequency that has a false alarm probability < 0.1 and power > highest fap
    faps = list(periodogram.false_alarm_probability(fap_levels))
    try: 
        best_frequency = freq[np.argmax(power)]
        best_power = power[np.argmax(power)]
        T = 1 / (float(best_frequency))
        phased_dates = np.mod(mjds, T) / T  # Got this from feets package documentation
        phased_dates_cycle_2 = phased_dates + 1
        
        if best_power > 0.1 and best_power > faps[(len(faps)-1)]:
            return best_power

    except ValueError:
        return 'NaN'
        pass


def test_mod_ls(data_file):
    #Import(s)
    import numpy as np
    import pandas as pd
    
    #Action
    df = pd.read_csv(data_file)
    mjds = np.array(df['mjd'])
    mags = np.array(df['mag'])
    magerrs = np.array(df['magerr'])
    
    return updated_ls(mjds=mjds,mags=mags,magerrs=magerrs,min_period=2,fap_levels=[0.001],out_type='')

print(test_mod_ls(data_file='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/2MASS_J20464241+4345102_r.csv'))
