def q_timeline(ids, data_temp, qs, plot_path):
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from matplotlib.patches import ConnectionPatch
    from matplotlib.ticker import AutoMinorLocator
    from astropy.timeseries import LombScargle
    
    # Action
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    fig = plt.figure(figsize=(9.5, 4))
    specs = fig.add_gridspec(2, 3)

    # Plot numberline
    n_ax = fig.add_subplot(specs[1, :])
    n_ax.spines['right'].set_color('none')
    n_ax.spines['left'].set_color('none')
    n_ax.yaxis.set_major_locator(ticker.NullLocator())
    n_ax.spines['top'].set_color('none')
    n_ax.patch.set_alpha(0.0)
    n_ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    n_ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    n_ax.set_xlim(0,1)
    n_ax.set_xlabel('Q')

    ### ADJUST SOME THINGS
    first_ax= fig.add_subplot(specs[0,0])
    second_ax = fig.add_subplot(specs[0,1])
    third_ax = fig.add_subplot(specs[0,2])
    for i in [first_ax, second_ax, third_ax]:
        i.legend(loc='upper right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
        i.tick_params(axis='both',which='both',direction='in')
        i.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        i.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        i.xaxis.set_minor_locator(AutoMinorLocator())
        i.yaxis.set_minor_locator(AutoMinorLocator())
        i.set_ylabel('r (mag)')
        i.invert_yaxis()

    first_ax.set_xlabel('Phase')
    second_ax.set_xlabel('Date (MJD)')
    third_ax.set_xlabel('Date (MJD)')


    ### Connected arrows

    def connect_to_q(q, xB):
        xyA = (q, 0)
        xyB = (xB, 0.8)
        coordsA = 'data'
        coordsB = 'data'
        con = ConnectionPatch(xyA, xyB, coordsA, coordsB,
                        arrowstyle="<|-|>",
                        mutation_scale=5, fc="black")
        n_ax.add_artist(con)

    connect_to_q(0.155, 0.14)
    connect_to_q(0.473, 0.5)
    connect_to_q(0.956, 0.856)
    
    ### GET AND PLOT DATA
    def get_data(id, data_temp):
        # Import(s)
        import pandas as pd 
        
        # Action
        data_file = data_temp.replace('+++', id)
        df = pd.read_csv(data_file)
        mjds = np.array(df['mjd'])
        mags = np.array(df['mag'])
        
        return mjds, mags
        
    ### FOLDED FIRST PLOT
    first_data = get_data('GDR1_2162251394832460672', data_temp)
    first_mjds = first_data[0]
    first_mags = first_data[1]
    
    periodogram = LombScargle(first_mjds,first_mags)
    all_freqs, all_powers = periodogram.autopower(method='fast',minimum_frequency=1/250,maximum_frequency=1/0.525)
    all_periods = 1/all_freqs
    best_per = all_periods[np.argmax(all_powers)]
    
    phased_dates = (np.mod(first_mjds, best_per))/best_per 
    phased_dates_cycle_2 = phased_dates + 1
    
    first_ax.scatter(phased_dates, first_mags, s=2, color = "#408ee0")
    first_ax.scatter(phased_dates_cycle_2, first_mags, s=2, color = "#408ee0")
    
    ## THE REST OF THE PLOTS
    second_data = get_data('2163140315626343168', data_temp)
    second_mjds = second_data[0]
    second_mags = second_data[1]
    second_ax.scatter(second_mjds, second_mags, color='#408ee0', edgecolor='black', linewidths=0.1)    
    
    third_data = get_data('FHK_348', data_temp)
    third_mjds = third_data[0]
    third_mags = third_data[1]
    third_ax.scatter(third_mjds, second_mags, color='#408ee0', edgecolor='black', linewidths=0.1)    

    plt.subplots_adjust(hspace=0.05, wspace=0.35)

    plt.show()
    #plt.savefig(plot_path, bbox_inches='tight', dpi=700)
    
q_timeline([0.12, 0.47, 0.93])