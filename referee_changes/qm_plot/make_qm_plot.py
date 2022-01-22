from numpy.core.fromnumeric import reshape


def make_qm_plot(results_file, q_bounds, plot_path): 
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    import seaborn as sns
    from sklearn.preprocessing import MinMaxScaler
    
    # Action
    results_df = pd.read_csv(results_file)
    qs = np.array(results_df['Q'])
    ms = np.array(results_df['M'])
    nus = np.array(results_df['NU'])
    #scaled_nus = 1000*(np.array(results_df['NU']))**1.25

    scaler = MinMaxScaler(feature_range=(10,60))
    reshaped_nus = nus.reshape(-1,1)
    scaler.fit(reshaped_nus)
    scaled_nus = scaler.transform(reshaped_nus)

    classes = np.array(results_df['primary_class'])
    
    # Split into classes
    def return_arrays(class_name):
        arrays_arr = [qs, ms, scaled_nus]
        return [i[np.where(classes==class_name)] for i in arrays_arr]
    
    p_qs, p_ms, p_nus = return_arrays('p')
    mp_qs, mp_ms, mp_nus = return_arrays('mp')
    qpd_qs, qpd_ms, qpd_nus = return_arrays('qpd')
    apd_qs, apd_ms, apd_nus = return_arrays('apd')
    qps_qs, qps_ms, qps_nus = return_arrays('qps')
    b_qs, b_ms, b_nus = return_arrays('b')
    l_qs, l_ms, l_nus = return_arrays('l')
    s_qs, s_ms, s_nus = return_arrays('s')
    
    # Make plot
    
    # Initial settings
    plt.rcParams['font.family'] = 'serif'
    fig, ax = plt.subplots(figsize=(7,5))
    
    # Plot classes
    ax.scatter(p_qs,p_ms,marker='s',color='C0', s=p_nus,label='Periodic',edgecolor='black',linewidths=0.5)
    
    # Quasi sym 
    ax.scatter(qps_qs,qps_ms,marker='o', color='C1', s=qps_nus,label='Quasi-periodic',edgecolor='black',linewidths=0.5)

    # Stochastic
    ax.scatter(s_qs,s_ms,marker='*', color='C2', s=s_nus,label='Stochastic',edgecolor='black',linewidths=0.5)
    
    # Aper dipper
    ax.scatter(apd_qs,apd_ms,marker='D', color='C3', s=apd_nus,label='Aperiodic dipper',edgecolor='black',linewidths=0.5)

    # Quasi periodic dippers
    ax.scatter(qpd_qs,qpd_ms,marker='h', color='C4', s=qpd_nus,label='Quasi-periodic dipper',edgecolor='black',linewidths=0.5)

    # Bursters
    ax.scatter(b_qs,b_ms,marker='d', color='C5', s=b_nus,label='Burster',edgecolor='black',linewidths=0.5)
    
    # Long timescales
    ax.scatter(l_qs,l_ms,marker='p', color='black', s=l_nus,label='Long-timescale', edgecolor='black',linewidths=0.5)
    
    # Multi periodics
    ax.scatter(mp_qs, mp_ms, marker='P', color='C6', s=mp_nus, label='Multi-periodic', edgecolor='black',linewidths=0.5)
    # Touch up settings
    
    plt.figtext(0.295,0.9,'Periodic',ha='center',fontsize='small')
    plt.figtext(0.64,0.9,'Quasi-Periodic',ha='center',fontsize='small')
    plt.figtext(0.84,0.9,'Aperiodic',ha='center',fontsize='small')

    plt.figtext(0.91,0.685,'Bursting',rotation=270,fontsize='small')
    plt.figtext(0.91,0.448,'Symmetric',rotation=270,fontsize='small')
    plt.figtext(0.91,0.24,'Dipping',rotation=270,fontsize='small')
    
    ax.axhline(y=0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    ax.axhline(y=-0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    ax.axvline(x=q_bounds[0],ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    ax.axvline(x=q_bounds[1],ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    ax.set(xlim=(0,1), ylim=(-1,1.1), xlabel='Quasi-Periodicity (Q)', ylabel='Flux Asymmetry (M)')
    
    ax.invert_yaxis()
   
    
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.tick_params(axis='both',which='both',direction='in')
    ax.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    ax.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)

    
    lgnd = plt.legend(loc='lower left',fontsize=6,fancybox=False,edgecolor='black',shadow=False,ncol=2)
    
    for i in range(8): 
        lgnd.legendHandles[i]._sizes = [18]
    
    #plt.show()
    plt.savefig('updated_qm_plot.png', dpi=300, bbox_inches='tight')

def make_scaled_qm_plot(results_file, q_bounds, plot_path): 
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    import seaborn as sns
    from sklearn.preprocessing import MinMaxScaler
    
    # Action
    results_df = pd.read_csv(results_file)
    qs = np.array(results_df['Q'])
    ms = np.array(results_df['M'])
    nus = np.array(results_df['NU'])
    #scaled_nus = 1000*(np.array(results_df['NU']))**1.25

    simple_scaler = MinMaxScaler(feature_range=(0,1))
    reshaped_qs = qs.reshape(-1,1)
    simple_scaler.fit(reshaped_qs)
    scaled_qs = simple_scaler.transform(reshaped_qs)

    scaled_q_bounds = simple_scaler.transform(np.reshape(q_bounds, (-1,1))).flatten()
    
    q_bounds = scaled_q_bounds
    qs = scaled_qs

    scaler = MinMaxScaler(feature_range=(10,60))
    reshaped_nus = nus.reshape(-1,1)
    scaler.fit(reshaped_nus)
    scaled_nus = scaler.transform(reshaped_nus)

    classes = np.array(results_df['primary_class'])
    
    # Split into classes
    def return_arrays(class_name):
        arrays_arr = [scaled_qs, ms, scaled_nus]
        return [i[np.where(classes==class_name)] for i in arrays_arr]
    
    p_qs, p_ms, p_nus = return_arrays('p')
    mp_qs, mp_ms, mp_nus = return_arrays('mp')
    qpd_qs, qpd_ms, qpd_nus = return_arrays('qpd')
    apd_qs, apd_ms, apd_nus = return_arrays('apd')
    qps_qs, qps_ms, qps_nus = return_arrays('qps')
    b_qs, b_ms, b_nus = return_arrays('b')
    l_qs, l_ms, l_nus = return_arrays('l')
    s_qs, s_ms, s_nus = return_arrays('s')
    
    # Make plot
    
    # Initial settings
    plt.rcParams['font.family'] = 'serif'
    fig, ax = plt.subplots(figsize=(7,5))
    
    # Plot classes
    ax.scatter(p_qs,p_ms,marker='s',color='C0', s=p_nus,label='Periodic',edgecolor='black',linewidths=0.5)
    
    # Quasi sym 
    ax.scatter(qps_qs,qps_ms,marker='o', color='C1', s=qps_nus,label='Quasi-periodic',edgecolor='black',linewidths=0.5)

    # Stochastic
    ax.scatter(s_qs,s_ms,marker='*', color='C2', s=s_nus,label='Stochastic',edgecolor='black',linewidths=0.5)
    
    # Aper dipper
    ax.scatter(apd_qs,apd_ms,marker='D', color='C3', s=apd_nus,label='Aperiodic dipper',edgecolor='black',linewidths=0.5)

    # Quasi periodic dippers
    ax.scatter(qpd_qs,qpd_ms,marker='h', color='C4', s=qpd_nus,label='Quasi-periodic dipper',edgecolor='black',linewidths=0.5)

    # Bursters
    ax.scatter(b_qs,b_ms,marker='d', color='C5', s=b_nus,label='Burster',edgecolor='black',linewidths=0.5)
    
    # Long timescales
    ax.scatter(l_qs,l_ms,marker='p', color='black', s=l_nus,label='Long-timescale', edgecolor='black',linewidths=0.5)
    
    # Multi periodics
    ax.scatter(mp_qs, mp_ms, marker='P', color='C6', s=mp_nus, label='Multi-periodic', edgecolor='black',linewidths=0.5)
    # Touch up settings
    
    plt.figtext(0.295,0.9,'Periodic',ha='center',fontsize='small')
    plt.figtext(0.64,0.9,'Quasi-Periodic',ha='center',fontsize='small')
    plt.figtext(0.84,0.9,'Aperiodic',ha='center',fontsize='small')

    plt.figtext(0.91,0.685,'Bursting',rotation=270,fontsize='small')
    plt.figtext(0.91,0.448,'Symmetric',rotation=270,fontsize='small')
    plt.figtext(0.91,0.24,'Dipping',rotation=270,fontsize='small')
    
    ax.axhline(y=0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    ax.axhline(y=-0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    ax.axvline(x=q_bounds[0],ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    ax.axvline(x=q_bounds[1],ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    ax.set(xlim=(0,1), ylim=(-1,1.1), xlabel='Quasi-Periodicity (Q)', ylabel='Flux Asymmetry (M)')
    
    ax.invert_yaxis()
   
    
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.tick_params(axis='both',which='both',direction='in')
    ax.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    ax.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)

    
    lgnd = plt.legend(loc='lower left',fontsize=6,fancybox=False,edgecolor='black',shadow=False,ncol=2)
    
    for i in range(8): 
        lgnd.legendHandles[i]._sizes = [18]
    
    #plt.show()
    plt.savefig('referee_changes/qm_plot/updated_qm_plot.png', dpi=300, bbox_inches='tight')
    

def change_scaling(results_file): 
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    import seaborn as sns
    from sklearn.preprocessing import MinMaxScaler
    
    # Action
    results_df = pd.read_csv(results_file)
    qs = np.array(results_df['Q'])
    ms = np.array(results_df['M'])
    nus = np.array(results_df['NU'])
    scaled_nus = 1000*(nus)**1.25
    log_nus = np.abs(np.log10(nus))

    fig, axs = plt.subplots(1, 4, figsize=(12,3))
    axs[0].hist(log_nus)
    axs[0].set(xlabel='new nu dist', ylabel='Frequency')

    better_log_nus = nus
    axs[1].hist(better_log_nus)
    axs[1].set(xlabel='better nu dist', ylabel='Frequency')

    axs[2].hist(scaled_nus)
    axs[2].set(xlabel='old scaled nu', ylabel='Frequency')
    
    scaler = MinMaxScaler(feature_range=(15,60))
    reshaped_nus = better_log_nus.reshape(-1,1)
    scaler.fit(reshaped_nus)
    new_nus = scaler.transform(reshaped_nus)
    axs[3].hist(new_nus)
    axs[3].set(xlabel='new scaled nu', ylabel='Frequency')

    plt.show()
    #plt.savefig('referee_changes/qm_plot/updated_qm_plot.png', dpi=250)


    '''
    plt.scatter(qs, ms, marker='o',color='C0', s=scaled_nus,label='All Classes',edgecolor='black',linewidths=0.5)

    plt.figtext(0.295,0.9,'Periodic',ha='center',fontsize='small')
    plt.figtext(0.64,0.9,'Quasi-Periodic',ha='center',fontsize='small')
    plt.figtext(0.848,0.9,'Aperiodic',ha='center',fontsize='small')

    plt.figtext(0.91,0.685,'Bursting',rotation=270,fontsize='small')
    plt.figtext(0.91,0.448,'Symmetric',rotation=270,fontsize='small')
    plt.figtext(0.91,0.24,'Dipping',rotation=270,fontsize='small')

    plt.axhline(y=0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axhline(y=-0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=0.47,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=0.85,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.xlim(0,1)
    plt.ylim(-1,1.1)
    plt.gca().invert_yaxis()
    plt.xlabel('Quasi-Periodicity (Q)')
    plt.ylabel('Flux Asymmetry (M)')
    plt.axes().yaxis.set_minor_locator(MultipleLocator(0.05))
    plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
    plt.tick_params(axis='both',which='both',direction='in')
    plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)

    lgnd = plt.legend(loc='lower left',fontsize=6,fancybox=False,edgecolor='black',shadow=False,ncol=2)

    plt.show()
    '''


results_file = './referee_changes/JANUARY_FINAL_RESULTS/JANUARY_FINAL_RESULTS.csv'

#change_scaling(results_file)
#make_qm_plot(results_file, [0.45,0.87], '')
make_scaled_qm_plot(results_file, [0.45,0.87], '')