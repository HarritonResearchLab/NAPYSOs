"""Functions for calculating Quasi-Periodicity Q and Flux Asymmetry M from
Bredall+ 2020, adapted from Cody+ 2014"""
import numpy as np
from astropy.convolution import Box1DKernel, convolve
from astropy.timeseries import LombScargle
from pandas import Series, concat, date_range, read_csv, to_datetime
from scipy.interpolate import interp1d


def best_per(JD: np.array, mag: np.array, err: np.array) -> float:
    r"""Calculate frequency with the highest power via Astropy LombScargle.
    Parameters
    ----------
    JD : np.array
        The julian dates for the lightcurve
    mag : np.array
        The magnitude measurements for the lightcurve
    err : np.array
    Returns
        The errors/uncertainties for the magnitude measurements
    -------
    out : float
        The best period, in days, found after filtering out the window function
        among other cuts...
    Notes
    -----
    Because ASAS-SN is a ground-based telescope network, the cadence alone can
    often result in frequencies with very high powers when fed into a standard
    Astropy LombScargle function. As such, we create a window function
    (Dirac Comb) out of the original light curve, i.e. values of 1 when an
    observation takes place and values of 0 elsewhere. We utilize Astropy
    LombScargle to find the periodograms of both the original curve and the
    Dirac Comb window function. We find all frequencies in the window function
    with a power greater than two standard deviations above the mean and ignore
    any peaks in the periodogram of the light curve within 0.03 Hz of these
    peaks in the window function, as well as any peaks with a false-alarm
    probability greater than 0.1.
    """
    # Generate a Dirac Comb, our window function
    time = to_datetime(JD, unit="D", origin="julian")
    time_not_obs = date_range(time.min(), time.max(), periods=1000)
    base = Series(np.zeros(len(time_not_obs)), index=time_not_obs)
    teeth = Series(np.ones(len(time)), index=time)
    dirac_comb = concat([base, teeth]).sort_index()

    # We only want to look for signals that we have 3 periods observed.
    maxp = (JD.max() - JD.min()) / 3
    minf = 1/maxp

    # First, we find the periodigram of the window function
    JD_W = dirac_comb.index.to_julian_date()
    mag_W = dirac_comb.values
    periodogram_W = LombScargle(JD_W, mag_W)
    freq_W, power_W = periodogram_W.autopower(
        minimum_frequency=minf, maximum_frequency=1/0.1
    )

    # Now, for the original lightcurve
    periodogram = LombScargle(JD, mag)
    freq, power = periodogram.autopower(
        minimum_frequency=minf, maximum_frequency=1/0.1
    )

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
    idx = np.argmax(power)
    best_fap = faps[idx]
    while best_fap > 0.1:
        power = np.delete(power, idx)
        idx = np.argmax(power)
        best_fap = faps[idx]
    best_freq = freq[idx]

    # Periods are more human-friendly :)
    return 1 / best_freq


def quas_per(JD, mag, err):
    r"""Calculate quasi-periodicity Q for a given lightcurve.
    ----------
    JD : np.array
        The julian dates for the lightcurve
    mag : np.array
        The magnitude measurements for the lightcurve
    err : np.array
        The errors/uncertainties for the magnitude measurements
    Returns
    -------
    out : float
        The quasi-periodicity of the lightcurve, typically 0 < Q < 1,
        with Q->0 => periodic and Q->1 => stocastic.
    Notes
    -----
    Quasi-periodicity requires knowledge of the estimated photometric
    uncertainty. What we did was grab 115,000 lightcurves in ASAS-SN for the
    corresponding filter (i.e. g or V), calculated their mean magnitude and
    stdev, and stored that in our reference file. Then, as you can see below,
    we did a rolling mean of the standard deviations as a function of mean_mag,
    allowing us to essentially get the typical variation for a given magnitude
    in a given filter in ASAS-SN. You may have some other way to get the
    photometric uncertainty; if so, adjust the code as needed, and store
    your photometric uncertainty in the variable 'sig'.
    """
    # Find the best period for the lc
    per = best_per(JD, mag, err)

    '''SIG CALCULATION AREA'''
    # First, we need a dataframe consiting of random lightcurves sorted
    # by magnitudes and their corresponding standard deviation
    rand_path = "PATH/TO/REFERENCE/FILE.csv"
    rand_filt = read_csv(
        rand_path, sep="\t", float_precision="high"
    ).loc[:, ["mean_mag", "std"]].sort_values("mean_mag").dropna()

    # Now we find a rolling average of the standard devation as a function
    # of magnitude. This is a little confusing; when "mean_mag" is in quotes,
    # this is referring to the mean magnitude of the lightcurve. When
    # we use the "mean()" method, we are taking the average std at a given
    # grouping of 500 mean magnitudes.
    std_per_mag = rand_filt.rolling(500, min_periods=1, on="mean_mag").mean()
    interp_mag = std_per_mag.loc[:, "mean_mag"]
    interp_std = std_per_mag.loc[:, "std"]

    # We interpolate this to get a nice continuous std as a function of mag
    phot_err = interp1d(interp_mag, interp_std, fill_value="extrapolate")
    sig = phot_err(np.nanmean(mag))
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