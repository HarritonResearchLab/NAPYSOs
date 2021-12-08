def updated_qm_plot(non_periodics,periodics,qm_results,q_bounds):
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.ticker import MultipleLocator
    import seaborn as sns

    # Action
    
    qm_results_df = pd.read_csv(qm_results)
    all_ids = np.array(qm_results_df['ID'])
    all_qs = np.array(qm_results_df['Q'])
    all_ms = np.array(qm_results_df['M'])    
    all_nus = np.array(qm_results_df['NU'])
    
    non_periodics_df = pd.read_csv(non_periodics)
    non_periodics_file_names = np.array(non_periodics_df['filename'])
    non_periodics_classes = np.array(non_periodics_df['class'])
    non_periodics_ids = np.array([])
    non_periodic_qs = np.array([])
    non_periodic_ms = np.array([])
    non_periodic_nus = np.array([])
    
    for filename in non_periodics_file_names:
        non_periodic_id = filename.split(':')[2].replace('.png','')
        non_periodics_ids = np.append(non_periodics_ids,non_periodic_id)
        for all_id, q, m, nu in zip(all_ids,all_qs,all_ms,all_nus):
            if non_periodic_id==all_id: 
                non_periodic_qs = np.append(non_periodic_qs,q)
                non_periodic_ms = np.append(non_periodic_ms,m)
                non_periodic_nus = np.append(non_periodic_nus,nu)
                break
            
    periodic_df = pd.read_csv(periodics)
    periodic_q = np.array(periodic_df['Q'])
    periodic_m = np.array(periodic_df['M'])
    periodic_nu = np.array(periodic_df['NU'])
                
    ### MODIFY NUS FOR SCALING PURPOSES
    periodic_nu = 1000*(periodic_nu**(2.5/2))
    non_periodic_nus = 1000*(non_periodic_nus**(2.5/2))
    
    lower_q_bound = q_bounds[0]
    upper_q_bound = q_bounds[1]
    
    ### SORT OTHER CLASSES
    #non_periodics_classes
    
    # Periodic sym 

    # Quasi periodic 
    quasi_periodic_mask = np.where(non_periodics_classes=='qps')
    quasi_q = non_periodic_qs[quasi_periodic_mask]
    quasi_m = non_periodic_ms[quasi_periodic_mask]
    quasi_nu = non_periodic_nus[quasi_periodic_mask]

    # Stochastic 
    stochastic_mask = np.where(non_periodics_classes=='s')
    stoch_q = non_periodic_qs[stochastic_mask]
    stoch_m = non_periodic_ms[stochastic_mask]
    stoch_nu = non_periodic_nus[stochastic_mask]

    # Aperiodic dippers
    aper_dipper_mask = np.where(non_periodics_classes=='apd')
    aper_dip_q = non_periodic_qs[aper_dipper_mask]
    aper_dip_m = non_periodic_ms[aper_dipper_mask]
    aper_dip_nu = non_periodic_nus[aper_dipper_mask]

    # Quasi periodic dippers 
    quasi_per_dip_mask = np.where(non_periodics_classes=='qpd')
    quasi_dipper_q = non_periodic_qs[quasi_per_dip_mask]
    quasi_dipper_m = non_periodic_ms[quasi_per_dip_mask]
    quasi_dipper_nu = non_periodic_nus[quasi_per_dip_mask]

    # Bursters
    burster_mask = np.where(non_periodics_classes=='b')
    burster_q = non_periodic_qs[burster_mask]
    burster_m = non_periodic_ms[burster_mask]
    burster_nu = non_periodic_nus[burster_mask]
    
    # Long timescales
    l_mask = np.where(non_periodics_classes=='l')
    l_q = non_periodic_qs[l_mask]
    l_m = non_periodic_ms[l_mask]
    l_nu = non_periodic_nus[l_mask]
    
    # Plot
    sns.color_palette()
    plt.rcParams['font.family'] = 'serif'
    plt.gca().tick_params(axis='both',which='both',direction='in')

    # Plot classes 
    #plt.scatter(marker='o',color='',s=,label='')
    
    # Periodic symmetric 
    plt.scatter(periodic_q,periodic_m,marker='o',s=periodic_nu,label='Periodic',edgecolor='black',linewidths=0.5)

    # Quasi sym 
    plt.scatter(quasi_q,quasi_m,marker='o',s=quasi_nu,label='Quasi-periodic',edgecolor='black',linewidths=0.5)

    # Stochastic
    plt.scatter(stoch_q,stoch_m,marker='o',s=stoch_nu,label='Stochastic',edgecolor='black',linewidths=0.5)

    # Aper dipper
    plt.scatter(aper_dip_q,aper_dip_m,marker='o',s=aper_dip_nu,label='Aperiodic dipper',edgecolor='black',linewidths=0.5)

    # Quasi periodic dippers
    plt.scatter(quasi_dipper_q,quasi_dipper_m,marker='o',s=quasi_dipper_nu,label='Quasi-periodic dipper',edgecolor='black',linewidths=0.5)

    # Bursters
    plt.scatter(burster_q,burster_m,marker='o',s=burster_nu,label='Burster',edgecolor='black',linewidths=0.5)
    
    # Long timescales
    plt.scatter(l_q,l_m,marker='o',s=l_nu,label='Long-timescale',edgecolor='black',color='black',linewidths=0.5)
    
    # Resume with rest of plot

    plt.figtext(0.29,0.9,'Periodic',ha='center',fontsize='small')
    plt.figtext(0.6275,0.9,'Quasi-Periodic',ha='center',fontsize='small')
    plt.figtext(0.84,0.9,'Aperiodic',ha='center',fontsize='small')

    plt.figtext(0.91,0.685,'Bursting',rotation=270,fontsize='small')
    plt.figtext(0.91,0.43,'Symmetric',rotation=270,fontsize='small')
    plt.figtext(0.91,0.22,'Dipping',rotation=270,fontsize='small')

    plt.axhline(y=0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axhline(y=-0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=lower_q_bound,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=upper_q_bound,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
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
    for i in range(7): 
        lgnd.legendHandles[i]._sizes = [15]
    
    plt.show()
    #plt.savefig('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/qm_classes_alpha/qm_vs_alpha.png',dpi=500)
    
    


non_periodics = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/qm_classes_alpha/non_periodic_qm_classes.csv'
periodics = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/qm_classes_alpha/qm_periodic_classes.csv'
qm_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results.csv'

updated_qm_plot(non_periodics,periodics,qm_results,q_bounds=[0.45,0.85])


'/home/thaddaeus/Random/Playing with Astropy/Unfolded.png'