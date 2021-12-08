def regression_plot(x, y): 
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import linregress
    from matplotlib.ticker import AutoMinorLocator
    
    # Action
    
    x_mask = np.logical_and(x>=0, y<=1)
    x, y = x[x_mask], y[x_mask]
    y_mask = np.logical_and(y>=0, y<=1)
    x, y = x[y_mask], y[y_mask]
    
    ols = linregress(x, y)
    
    m = ols[0]
    b = ols[1]
    r_sq = ols[2]**2
    
    # Make plot
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.figsize'] = (5, 5)
    
    plt.scatter(x, y, color='#408ee0', edgecolors='black', linewidths=0.5)
    
    x_arr = np.array([x.min, x.max])
    
    ols_label = 'm: ' + str(round(m, 1)) + '\n' + r'$R^2$' + ': ' + str(round(r_sq, 2))
    
    plt.plot(x_arr, m*x_arr+b, color='black', linewidths=0.5, ls='--', label=ols_label)
    
    plt.axline(xy1=(0, 0), slope=1, color='grey', ls='--', label='1:1')
    
    #plt.legend()
    
    plt.tick_params(axis='both',which='both',direction='in')
    plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator())    
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
    
    plt.show()
    
    return m, b, r_sq
    
    
'/home/thaddaeus/FMU/HRL/LAH2.0/efforts/omicron/kepler_regs/data/Keppler Q.csv'


import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 10)

