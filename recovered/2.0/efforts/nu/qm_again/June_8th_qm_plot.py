def june8_qm_plot(results_file, q_bounds, plot_path): 
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    import seaborn as sns
    
    # Action
    results_df = pd.read_csv(results_file)
    qs = np.array(results_df['Q'])
    ms = np.array(results_df['M'])
    scaled_nus = 1000*(np.array(results_df['NU']))**1.25
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
    plt.gca().tick_params(axis='both',which='both',direction='in')
    
    # Plot classes
    plt.scatter(p_qs,p_ms,marker='o',color='C0', s=p_nus,label='Periodic',edgecolor='black',linewidths=0.5)

    # Quasi sym 
    plt.scatter(qps_qs,qps_ms,marker='o', color='C1', s=qps_nus,label='Quasi-periodic',edgecolor='black',linewidths=0.5)

    # Stochastic
    plt.scatter(s_qs,s_ms,marker='o', color='C2', s=s_nus,label='Stochastic',edgecolor='black',linewidths=0.5)

    # Aper dipper
    plt.scatter(apd_qs,apd_ms,marker='o', color='C3', s=apd_nus,label='Aperiodic dipper',edgecolor='black',linewidths=0.5)

    # Quasi periodic dippers
    plt.scatter(qpd_qs,qpd_ms,marker='o', color='C4', s=qpd_nus,label='Quasi-periodic dipper',edgecolor='black',linewidths=0.5)

    # Bursters
    plt.scatter(b_qs,b_ms,marker='o', color='C5', s=b_nus,label='Burster',edgecolor='black',linewidths=0.5)
    
    # Long timescales
    plt.scatter(l_qs,l_ms,marker='o', color='black', s=l_nus,label='Long-timescale', edgecolor='black',linewidths=0.5)
    
    # Multi periodics
    plt.scatter(mp_qs, mp_ms, marker='o', color='C6', s=mp_nus, label='Multi-periodic', edgecolor='black',linewidths=0.5)
    # Touch up settings
    
    plt.figtext(0.295,0.9,'Periodic',ha='center',fontsize='small')
    plt.figtext(0.64,0.9,'Quasi-Periodic',ha='center',fontsize='small')
    plt.figtext(0.848,0.9,'Aperiodic',ha='center',fontsize='small')

    plt.figtext(0.91,0.685,'Bursting',rotation=270,fontsize='small')
    plt.figtext(0.91,0.448,'Symmetric',rotation=270,fontsize='small')
    plt.figtext(0.91,0.24,'Dipping',rotation=270,fontsize='small')

    plt.axhline(y=0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axhline(y=-0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=q_bounds[0],ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=q_bounds[1],ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
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
    
    for i in range(8): 
        lgnd.legendHandles[i]._sizes = [18]
        
    #plt.show()
    plt.savefig(plot_path,dpi=300, bbox_inches='tight')
    
results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/september21/SEPTEMBER_FINAL_RESULTS/merged_results.csv'
q_bounds = [0.45, 0.87]
plot_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/september21/20-26/plots_for_paper/qm_plot.png'

june8_qm_plot(results_file, q_bounds, plot_path)