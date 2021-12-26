def plot_rotation_periods(ref_file):
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from scipy.optimize import curve_fit
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator


    # Action

    df = pd.read_csv(ref_file)
    periods = list(df['PER1'])
    powers = list(df['POW1'])
    faps = list(df['99%_FAP'])
    qs = list(df['Q'])

    good_periods = []
    for period, pow, fap, q in zip(periods,powers,faps,qs):
        if pow > fap: 
            if q < 0.45:
                if period<18: 
                    good_periods.append(period)

    periods = np.array(periods)
    
    # Fitting
    def gauss(x,mu,sigma,A):
        return A*np.exp(-(x-mu)**2/2/sigma**2)

    def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
        return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)


    # Plot
    plt.rcParams['font.family'] = 'serif'
    plt.hist(good_periods,bins=20,color='#408ee0',edgecolor='black',)
    plt.tick_params(axis='both',which='both',direction='in')
    plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
    plt.gca().yaxis.set_minor_locator(MultipleLocator(1))
    plt.xlabel('Period (d)')
    plt.ylabel('Frequency')
    plt.show()


# ACTUALLY RUN IT
plot_rotation_periods('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results.csv')