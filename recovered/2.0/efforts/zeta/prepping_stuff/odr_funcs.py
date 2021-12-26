

def fit_twocolor_odr(band1, band2, band1_err, band2_err,  n_bootstrap = 1000, xyswitch = False, p_guess = None):
    '''Fits a straight line to a single CMD, using a weighted orthogonal least squares algorithm (ODR).
    
    Parameters
    ----------
    data1 : np.array
        single lightcurve of band 1 in magnitudes
    data2 : np.array
        single lightcurve of band 2 in magnitudes      
    data1_error : np.array
        error on data points of band 1 in magnitudes    
    data2_error : np.array
        error on data points of band 2 in magnitudes
    dataset : np.ndarray
        data collection for one detected source
    index : integer
        the index of the dataset within the data structure
    p_guess : tuple
        initial fit parameters derived from fit_twocolor
    outroot : string or None
        dictionary where to save the plot, set to `None` for no plotting
    n_bootstrap : integer or None
        how many bootstrap trials, set to `None` for no bootstrapping
    xyswitch : boolean
        if the X and Y axis will be switched for the fit or not. This has nothing to do with bisector fitting! The fitting algorithm used here takes care of errors in x and y simultaneously; the xyswitch is only for taking care of pathological cases where a vertical fitted line would occur without coordinate switching.
    
    Returns
    -------
    result : tuple
        contains output = fit parameters, bootstrap_output = results from the bootstrap, bootstrap_raw = the actual bootstrapped data, alpha = the fitted slope angle, sd_alpha = the error on the fitted slope angle, x_spread = the spread of the data along the fitted line (0.5*(90th percentile - 10th percentile)))
    '''
    #Imports
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    import scipy
    from scipy import odr, stats

    #Action
    #define the fitting function (in this case a straight line)
    def fitfunc(p, x):
        return p[0]*x + p[1]

    if p_guess is None:
        simple_cmd = cmd_slope_simple(band1, band2, band1_err, band2_err)
        if type(simple_cmd) != float:
            p_guess = list(simple_cmd)[0:2]
            if np.isfinite(p_guess[0])==True: # pathological case
                p_guess[0] = 0
            if np.isfinite(p_guess[1])==True: # pathological case
                p_guess[1] = np.mean(band1-band2)
    
            # define what the x and y data is:
            x_data = band1 - band2
            y_data = band1
            x_error = np.sqrt( band1_err**2 + band2_err**2 )
            y_error = band1_err
            if xyswitch:
                y_data, x_data = (x_data, y_data)
                y_error, x_error = (x_error, y_error)
            
            # load data into ODR
            data = scipy.odr.RealData(x=x_data, y=y_data, sx=x_error, sy=y_error)
            # tell ODR what the fitting function is:
            model = scipy.odr.Model(fitfunc)
            # now do the fit:
            fit = scipy.odr.ODR(data, model, p_guess, maxit=1000) 
            output = fit.run()
            
            p = output.beta # the fitted function parameters
            delta = output.delta # array of estimated errors in input variables
            eps   = output.eps # array of estimated errors in response variables
            #print output.stopreason[0]
            bootstrap_output = np.array([np.NaN, np.NaN, np.NaN, np.NaN])
            bootstrap_raw = (np.NaN, np.NaN, np.NaN)
            # calculate slope angle. This is vs. horizontal axis.
            alpha = math.atan(output.beta[0])
            # calculate error on slope angle by taking the mean difference of the angles derived from m+m_error and m-m_error.
            alpha_plus  = math.asin((output.beta[0]+output.sd_beta[0])/np.sqrt((output.beta[0]+output.sd_beta[0])**2 + 1**2))
            alpha_minus = math.asin((output.beta[0]-output.sd_beta[0])/np.sqrt((output.beta[0]-output.sd_beta[0])**2 + 1**2))
            sd_alpha = 0.5*( np.abs(alpha - alpha_plus) + np.abs(alpha - alpha_minus) ) 
            # define the spread along the fitted line. Use 90th and 10th quantile.
            # output.xplus and output.y are the x and y values of the projection of the original data onto the fit.
            # okay, first transform coordinate system so that x axis is along fit. To do this, first shift everything by -p[1] (this is -b), then rotate by -alpha. New x and y coordinates are then:
            #
            # |x'|   |cos(-alpha) -sin(-alpha)| | x |
            # |  | = |                        | |   |
            # |y'|   |sin(-alpha)  cos(-alpha)| |y-b|
            #
            x_new = math.cos(-alpha) * output.xplus - math.sin(-alpha)*(output.y - p[1])
            y_new = math.sin(-alpha) * output.xplus + math.cos(-alpha)*(output.y - p[1])
            # The y_new values are now essentially zero. (As they should.)
            # Now sort x_new and get 90th and 10th quantile:
            x_new.sort()
            x_spread = scipy.stats.mstats.mquantiles(x_new, prob=0.9)[0] - scipy.stats.mstats.mquantiles(x_new, prob=0.1)[0]
            #print x_spread
            
            if n_bootstrap is not None:
                #print('bootstrapping...')
                # take a random half of the data and do the fit (choosing without replacement, standard bootstrap). Do this a lot of times and construct a cumulative distribution function for the slope and the intercept of the fitted line.
                # now what I actually want is the slope angle a, not m.
                m = np.array([])
                b = np.array([])
                for i in np.arange(0, n_bootstrap):
                    indices = np.arange(0,len(x_data))
                    np.random.shuffle(indices)
                    ind = indices[0:int(round(len(x_data)/2,0))] # dividing by integer on purpose.
                    dat = scipy.odr.RealData(x=x_data[ind], y=y_data[ind], sx=x_error[ind], sy=y_error[ind])
                    fit = scipy.odr.ODR(dat, model, p_guess, maxit=5000,job=10) 
                    out = fit.run()
                    m = np.append(m, out.beta[0])
                    b = np.append(b, out.beta[1])
                
                a = np.arctan(m) # in radian
                # get median and symmetric 68% interval for m, b and alpha:
                m_median = np.median(m)
                m_down = np.sort(m)[ int(round(0.16*len(m))) ]
                m_up   = np.sort(m)[ int(round(0.84*len(m))) ]
                m_error = np.mean([abs(m_down-m_median), abs(m_up-m_median)])
                #print (m_median, m_up, m_down, m_error)
                b_median = np.median(b)
                b_down = np.sort(b)[ int(round(0.16*len(b))) ]
                b_up   = np.sort(b)[ int(round(0.84*len(b))) ]
                b_error = np.mean([abs(b_down-b_median), abs(b_up-b_median)])
                #print (b_median, b_up, b_down, b_error)
                a_median = np.median(a)
                a_down = np.sort(a)[ int(round(0.16*len(a))) ]
                a_up   = np.sort(a)[ int(round(0.84*len(a))) ]
                a_error = np.mean([abs(a_down-a_median), abs(a_up-a_median)])
                #print (b_median, b_up, b_down, b_error)
                
                #bootstrap_output = np.array([m_median, m_error, b_median, b_error, a_median, a_error])
                #bootstrap_raw = (m, b, a)
            
            #result = (output, bootstrap_output, bootstrap_raw, alpha, sd_alpha, x_spread)
            #return result
            
            slope_angle = math.degrees(a_median)
            slope_angle_error = math.degrees(a_error)

            return [slope_angle,slope_angle_error,x_spread]

def test_it(ID,csv_temp): 
    #Import(s)
    import pandas as pd
    import numpy as np

    #Action

    csv_file = csv_temp.replace('+++',ID)
    df = pd.read_csv(csv_file)

    g_mags = np.array(df['g_mag'])
    g_magerrs = np.array(df['g_magerr'])
    r_mags = np.array(df['r_mag'])
    r_magerrs = np.array(df['r_magerr'])

    """
    This is probably where you want to restrict the cmd data. 
    """

    paup_odr_results = fit_twocolor_odr(band1=g_mags,band2=r_mags,band1_err=g_magerrs,band2_err=r_magerrs)
    if paup_odr_results != None:
        return paup_odr_results
    else:
        return 'NaN'

def updated_cmd_routine(key):
    #Imports
    import pandas as pd
    import math
    from progress.bar import Bar
    
    #Action
    in_df = pd.read_csv(key)
    ids = list(in_df['ID']) 

    good_ids = []
    slope_angles = []
    slopes = []
    slope_errors = []
    spreads = []

    bar = Bar('Processing...',max=len(ids))

    for i ,id in enumerate(ids):    
        odr_results = test_it(ID=id,csv_temp='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/cmd_v2/+++.csv')
        if odr_results != 'NaN':
            slope_angle = odr_results[0]
            slope = math.tan(math.radians(slope_angle))
            slope_error = odr_results[1]

            spread = odr_results[2]
        
            good_ids.append(id)
            slope_angles.append(slope_angle)
            slopes.append(slope)
            slope_errors.append(slope_error)
            spreads.append(spread)

        bar.next()
    
    zipped_list = list(zip(good_ids,slopes,slope_angles,slope_errors,spreads))
    out_df = pd.DataFrame(data=zipped_list,columns=['ID','SLOPE','SLOPE_ANGLE','SLOPE_ANGLE_ERROR','80%_SPREAD'])
    out_df.to_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/epsilon/updated_odr_results/updated_odr_results.csv',index=False)
    bar.finish()
    
updated_cmd_routine(key='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/key.csv')