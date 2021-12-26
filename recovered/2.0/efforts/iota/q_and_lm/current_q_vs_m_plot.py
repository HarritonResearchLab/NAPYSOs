def fig_plot_2(data_file,q_bounds):
    '''
    Create figure two plot (from Bredall et al. 2020) from data file which has all the IDs, Qs, and Ms 
    '''
    
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.ticker import MultipleLocator
    import seaborn as sns

    # Action

    df = pd.read_csv(data_file)
    ids = list(df['ID'])
    q_raw = np.array(df['Q'])

    q_mask = np.logical_and(q_raw>0.0,q_raw<1.0)
    
    m = np.array(df['M'])[q_mask]
    q_raw = q_raw[q_mask]
    nus = np.array(df['NU'])[q_mask]
    nus = 1000*(nus**(2.5/2))
    
    # Transform Q
    '''
    q = (q_raw-min(q_raw))*(1/(max(q_raw)-min(q_raw)))
    lower_q_bound = (q_bounds[0]-min(q_raw))*(1/(max(q_raw)-min(q_raw)))
    upper_q_bound = (q_bounds[1]-min(q_raw))*(1/(max(q_raw)-min(q_raw)))
    '''
    lower_q_bound = q_bounds[0]
    upper_q_bound = q_bounds[1]
    q = q_raw

    # Get morph classes 

    # Six classes 

    sym_ms_mask = np.where(np.logical_and(0.25>m,m<0.25))
    quasi_q_mask = np.where(np.logical_and(q<0.85,q>0.45))
    dipper_mask = np.where(m>0.25)

    # Periodic sym 
    periodic_qs_mask = np.where(q<0.45)
    periodic_mask = np.intersect1d(periodic_qs_mask,sym_ms_mask)
    periodic_q = q[periodic_mask]
    periodic_m = m[periodic_mask]
    periodic_nu = nus[periodic_mask]

    # Quasi periodic 
    quasi_periodic_mask = np.intersect1d(quasi_q_mask,sym_ms_mask)
    quasi_q = q[quasi_periodic_mask]
    quasi_m = m[quasi_periodic_mask]
    quasi_nu = nus[quasi_periodic_mask]

    # Stochastic 
    stochastic_q_mask = np.where(q>0.85)
    stochastic_m_mask = np.where(np.logical_and(m>-0.25,m<0.25))
    stochastic_mask = np.intersect1d(stochastic_m_mask,stochastic_q_mask)
    stoch_q = q[stochastic_mask]
    stoch_m = m[stochastic_mask]
    stoch_nu = nus[stochastic_mask]

    # Aperiodic dippers
    aper_dipper_mask = np.intersect1d(stochastic_q_mask,dipper_mask)
    aper_dip_q = q[aper_dipper_mask]
    aper_dip_m = m[aper_dipper_mask]
    aper_dip_nu = nus[aper_dipper_mask]

    # Quasi periodic dippers 
    quasi_per_dip_mask = np.intersect1d(quasi_q_mask,dipper_mask)
    quasi_dipper_q = q[quasi_per_dip_mask]
    quasi_dipper_m = m[quasi_per_dip_mask]
    quasi_dipper_nu = nus[quasi_per_dip_mask]

    # Bursters
    burster_mask = np.where(m<-0.25)
    burster_q = q[burster_mask]
    burster_m = m[burster_mask]
    burster_nu = nus[burster_mask]

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

    plt.legend(loc='lower left',fontsize=6,fancybox=False,edgecolor='black',shadow=False,ncol=2)

    plt.show()

fig_plot_2(data_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results.csv',q_bounds=[0.45,0.85])