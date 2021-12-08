def fig_plot_2(data_file):
    '''
    Create figure two plot (from Bredall et al. 2020) from data file which has all the IDs, Qs, and Ms 
    '''
    
    # Import(s)
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.ticker import MultipleLocator

    # Action

    df = pd.read_csv(data_file)
    ids = list(df['ID'])
    q = list(df['Q'])
    m = list(df['M'])

    # Plot
    plt.rcParams['font.family'] = 'serif'

    plt.scatter(q,m,marker='.',color='slategray',s=2)

    plt.figtext(0.18,0.9,'Periodic',ha='center',fontsize='small')
    plt.figtext(0.5,0.9,'Quasi-Periodic',ha='center',fontsize='small')
    plt.figtext(0.845,0.9,'Aperiodic',ha='center',fontsize='small')

    plt.figtext(0.91,0.685,'Bursting',rotation=270,fontsize='small')
    plt.figtext(0.91,0.43,'Symmetric',rotation=270,fontsize='small')
    plt.figtext(0.91,0.22,'Dipping',rotation=270,fontsize='small')

    plt.axhline(y=0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axhline(y=-0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=0.15,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=0.85,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.xlim(0,1)
    plt.ylim(-1,1)
    plt.gca().invert_yaxis()
    plt.xlabel('Quasi-Periodicity (Q)')
    plt.ylabel('Flux Asymmetry (M)')
    plt.axes().yaxis.set_minor_locator(MultipleLocator(0.05))
    plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
    plt.show()


