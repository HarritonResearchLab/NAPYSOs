def do_stuff(z):
    #Import(s)
    import numpy as np
    from progress.bar import Bar
    import pandas as pd
    
    #Action
    slopes = []
    median_gs = []
    slope_errors = []
    ids = []

    with open('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/alpha/slope_vs_med_g/odr_results.csv','r') as f:
        for index, line in enumerate(f):
            if index != 0 and index %2 == 0:
                line = line.replace('\n','')
                line_list = line.split(',')
                id = line_list[0]
                slope = float(line_list[1])
                sd_beta = float(line_list[2])
                slope_error = (100*((z*sd_beta)/slope))
                if 0 < slope_error < 10 and slope>0:
                    g_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++_g.csv'.replace('+++',id)
                    g_mags = np.array(pd.read_csv(g_file)['mag'])
                    g_median = np.median(g_mags)
                    slopes.append(slope)
                    median_gs.append(g_median)
                    slope_errors.append(slope_error)
                    ids.append(id)
    
    return slopes, median_gs, slope_errors, ids

def plot_some():
    #Import(s)
    import matplotlib.pyplot as plt

    #Action

    z1 = do_stuff(z=1)
    z2 = do_stuff(z=1.645)
    z3 = do_stuff(z=1.96)
    z4 = do_stuff(z=2.58)

    for slope in z1[0]:
        if slope > 2.5:
            slope_index = z1[0].index(slope)
            print(z1[3][slope_index])

    #Plot
    plt.rcParams['font.family'] = 'serif'
    fig, axs = plt.subplots(2,2)


    axs[0,0].scatter(z1[0],z1[1],color='cornflowerblue',marker='.')
    axs[0,0].set_xlabel('Slope (z=1)')
    axs[0,0].set_ylabel('Median g')
    
    axs[0,1].scatter(z2[0],z2[1],color='cornflowerblue',marker='.')
    axs[0,1].set_xlabel('Slope (z=1.645)')
    axs[0,1].set_ylabel('Median g')

    axs[1,0].scatter(z3[0],z3[1],color='cornflowerblue',marker='.')
    axs[1,0].set_xlabel('Slope (z=1.96)')
    axs[1,0].set_ylabel('Median g')

    axs[1,1].scatter(z4[0],z4[1],color='cornflowerblue',marker='.')
    axs[1,1].set_xlabel('Slope (z=2.58)')
    axs[1,1].set_ylabel('Median g')
    
    fig.suptitle('Error formula: '+r'$100\cdot\frac{z\cdot \mathrm{SE}}{\mathrm{slope}}$')
    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    plt.show()



         
plot_some()