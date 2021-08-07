def Q(dates, mags, magerrs, timescale, sig_factor): 
    '''
    Parameters: 
    
    dates: array of time values associated with observations
    mags: array of magnitude values associated with observations
    magerrs: array of errors associated with magnitudes in mags array
    timescale: float representative of an object's timescale of variability
    sig_factor: multiplied against mean magerr, should be set to 1. initially
    
    Returns: 
    q: quasi-periodicity value
    resid_mags: array of residual mags
    '''
    
    # Import(s)
    import numpy as np
    from astropy.convolution import Box1DKernel, convolve
    
    # Convert to arrays
    mjd = np.array(mjd)
    mag = np.array(mag)
    magerr = np.array(magerr)
        
    # Calculate sig
    sig = sig_factor*np.mean(magerr)

    # Create the residual curve
    phase = mjd % timescale
    mag = mag[np.argsort(phase)]

    # We use three periods and extract the middle to prevent edge effects
    three_periods = np.concatenate((mag, mag, mag))
    boxcar = Box1DKernel(len(mag) // 4)
    smooth_mag = convolve(three_periods, boxcar)
    smooth_mag = smooth_mag[np.size(mag):2*np.size(mag)]

    resid_mag = mag - smooth_mag

    q = float((np.nanstd(resid_mag) ** 2 - sig ** 2)/(np.nanstd(mag) ** 2 - sig ** 2))

    return q, resid_mag
    
    
def M(mags):
    '''
    Parameters: 

    mags: array of magnitude values 
    
    Returns: 
    
    m: flux-asymmetry value
    '''
    
    # Import(s)
    from scipy import stats 
    import numpy as np
    
    # Calculate m 
    m = (np.mean([stats.mstats.mquantiles(mag,prob=0.9),stats.mstats.mquantiles(mag,prob=0.1)])-np.median(mag))/np.sqrt(((mag-mag.mean())**2).sum()/len(mag))
    
    return m
