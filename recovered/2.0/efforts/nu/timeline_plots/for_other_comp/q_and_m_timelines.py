def q_timeline(ids, data_temp, qs, plot_path, new_ids):
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from matplotlib.patches import ConnectionPatch
    from matplotlib.ticker import AutoMinorLocator
    from astropy.timeseries import LombScargle
    from matplotlib.ticker import MultipleLocator
    
    # Action
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    fig = plt.figure(figsize=(10, 4))
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
    n_ax.xaxis.set_major_locator(MultipleLocator(0.2))
    n_ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    
    ### ADJUST SOME THINGS
    first_ax= fig.add_subplot(specs[0,0])
    second_ax = fig.add_subplot(specs[0,1])
    third_ax = fig.add_subplot(specs[0,2])
    for i in [first_ax, second_ax, third_ax]:
        #i.legend(loc='upper right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
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

    connect_to_q(qs[0], 0.14)
    connect_to_q(qs[1], 0.5)
    connect_to_q(qs[2], 0.856)
    
    
    
    ### GET AND PLOT DATA
    def get_data(id, data_temp):
        # Import(s)
        import pandas as pd 
        
        # Action
        data_file = data_temp.replace('+++', id)
        df = pd.read_csv(data_file)
        mjds = np.array(df['mjd'])
        mags = np.array(df['mag'])
        magerrs = np.array(df['magerr'])
        
        return mjds, mags, magerrs
        
    ### FOLDED FIRST PLOT
    first_data = get_data(ids[0], data_temp)
    first_mjds = first_data[0]
    first_mags = first_data[1]
    first_magerrs = first_data[2]
    
    periodogram = LombScargle(first_mjds,first_mags)
    all_freqs, all_powers = periodogram.autopower(method='fast',minimum_frequency=1/250,maximum_frequency=1/0.525)
    all_periods = 1/all_freqs
    best_per = all_periods[np.argmax(all_powers)]
    
    phased_dates = (np.mod(first_mjds, best_per))/best_per 
    phased_dates_cycle_2 = phased_dates + 1
    
    first_ax.errorbar(phased_dates, first_mags, yerr = first_magerrs, color='#408ee0', lw=0, elinewidth=0.25, alpha=0.4)
    first_ax.errorbar(phased_dates_cycle_2, first_mags, yerr = first_magerrs, color='#408ee0', lw=0, elinewidth=0.25, alpha=0.4)
    
    first_ax.scatter(phased_dates, first_mags, s=6, color = "#408ee0", edgecolor='black', linewidths=0.2, marker='o')
    first_ax.scatter(phased_dates_cycle_2, first_mags, s=6, color = "#408ee0", edgecolor='black', linewidths=0.2, marker='o')
    
    first_ax.set_title(new_ids[0].replace('_',' '), fontsize=9)
    
    ## THE REST OF THE PLOTS
    second_data = get_data(ids[1], data_temp)
    second_mjds = second_data[0]
    second_mags = second_data[1]
    second_magerrs = second_data[2]
    
    second_ax.errorbar(second_mjds, second_mags, yerr=second_magerrs, color='#408ee0', lw=0, elinewidth=0.25, alpha=0.4)
    second_ax.scatter(second_mjds, second_mags, s=6, color='#408ee0', edgecolor='black', linewidths=0.2, marker='o')    
    second_ax.set_title(new_ids[1].replace('_',' '), fontsize=9)
    
    third_data = get_data(ids[2], data_temp)
    third_mjds = third_data[0]
    third_mags = third_data[1]
    third_magerrs = third_data[2]
    
    third_ax.errorbar(third_mjds, third_mags, yerr=third_magerrs, color='#408ee0', lw=0, elinewidth=0.25, alpha=0.4)
    third_ax.scatter(third_mjds, third_mags, s=6, color='#408ee0', edgecolor='black', linewidths=0.2, marker='o')    
    third_ax.set_title(new_ids[2].replace('_',' '), fontsize=9)
    
    plt.subplots_adjust(hspace=0.05, wspace=0.4)

    #plt.show()
    plt.savefig(plot_path, bbox_inches='tight', dpi=400)

def m_timeline(ids, data_temp, ms, plot_path, new_ids): 
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from matplotlib.patches import ConnectionPatch
    from matplotlib.ticker import AutoMinorLocator
    from matplotlib.ticker import MultipleLocator
    
    # Action
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    fig = plt.figure(figsize=(10, 4))
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
    n_ax.set_xlim(-1.3,1.)
    n_ax.set_xlabel('M')
    n_ax.xaxis.set_major_locator(MultipleLocator(0.4))
    n_ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    
    ### ADJUST SOME THINGS
    first_ax = fig.add_subplot(specs[0,0])
    second_ax = fig.add_subplot(specs[0,1])
    third_ax = fig.add_subplot(specs[0,2])
    for i in [first_ax, second_ax, third_ax]:
        i.tick_params(axis='both',which='both',direction='in')
        i.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        i.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        i.xaxis.set_minor_locator(AutoMinorLocator())
        i.yaxis.set_minor_locator(AutoMinorLocator())
        i.set_ylabel('r (mag)')
        i.invert_yaxis()
        i.set_xlabel('Date (MJD)')

    ### Connected arrows

    def connect_to_m(m, xB):
        xyA = (m, 0)
        xyB = (xB, 0.8)
        coordsA = 'data'
        coordsB = 'data'
        con = ConnectionPatch(xyA, xyB, coordsA, coordsB,
                        arrowstyle="<|-|>",
                        mutation_scale=5, fc="black")
        n_ax.add_artist(con)

    connect_to_m(ms[0], -1.)
    connect_to_m(ms[1], -0.165)
    connect_to_m(ms[2], 0.685)
    
    ### GET AND PLOT DATA
    def get_data(id, data_temp):
        # Import(s)
        import pandas as pd 
        
        # Action
        data_file = data_temp.replace('+++', id)
        df = pd.read_csv(data_file)
        mjds = np.array(df['mjd'])
        mags = np.array(df['mag'])
        magerrs = np.array(df['magerr'])
        
        return mjds, mags, magerrs
    
    first_data = get_data(ids[0], data_temp)
    first_mjds = first_data[0]
    first_mags = first_data[1]
    first_magerrs = first_data[2]
    
    first_ax.errorbar(first_mjds, first_mags, yerr=first_magerrs, color='#408ee0', lw=0, elinewidth=0.25, alpha=0.4)
    first_ax.scatter(first_mjds, first_mags, s=6, color='#408ee0', edgecolor='black', linewidths=0.2, marker='o')    
    first_ax.set_title(new_ids[0].replace('_',' '), fontsize=9)
        
    second_data = get_data(ids[1], data_temp)
    second_mjds = second_data[0]
    second_mags = second_data[1]
    second_magerrs = second_data[2]
    
    second_ax.errorbar(second_mjds, second_mags, yerr=second_magerrs, color='#408ee0', lw=0, elinewidth=0.25, alpha=0.4)
    second_ax.scatter(second_mjds, second_mags, s=6, color='#408ee0', edgecolor='black', linewidths=0.2, marker='o')    
    second_ax.set_title(new_ids[1].replace('_',' '), fontsize=9)
    
    third_data = get_data(ids[2], data_temp)
    third_mjds = third_data[0]
    third_mags = third_data[1]
    third_magerrs = third_data[2]
    
    third_ax.errorbar(third_mjds, third_mags, yerr=third_magerrs, color='#408ee0', lw=0, elinewidth=0.25, alpha=0.4)
    third_ax.scatter(third_mjds, third_mags, s=6, color='#408ee0', edgecolor='black', linewidths=0.2, marker='o')    
    third_ax.set_title(new_ids[2].replace('_',' '), fontsize=9)
    
    plt.subplots_adjust(hspace=0.05, wspace=0.4)

    #plt.show()
    plt.savefig(plot_path, bbox_inches='tight', dpi=400)

    
q_ids = ['GDR1_2162251394832460672', 'GDR1_2163140315626343168', 'FHK_348'] 
new_q_ids = ['FHK_446','LkHA_149', 'FHK_348']

data_temp = './recovered/2.0/data/AUGUST_5th/light_curves/+++_r.csv'   
qs = [0.155, 0.473, 0.956]
q_plot_path = './referee_changes/new_names/new_figs/q_timeline.png'

q_timeline(q_ids, data_temp, qs, q_plot_path, new_q_ids)

m_ids = ['GDR1_2162872928138757248', '2MASS_J20523394+4429168', 'FHK_577']
new_m_ids = ['FHK_267', '2MASS_J20523394+4429168', 'FHK_577']
ms = [-1.1991309270389323, 0.08255161814737327, 0.8993104086930837]
m_plot_path = './referee_changes/new_names/new_figs/m_timeline.png'

m_timeline(m_ids, data_temp, ms, m_plot_path, new_m_ids)