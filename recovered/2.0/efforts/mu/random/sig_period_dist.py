def sig_histogram(data_file, save_path):
    # Import(s)
    import numpy as np
    import pandas as pd 
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    
    # Action
    data_df = pd.read_csv(data_file)
    periods = np.array(data_df['PER1'])[np.where(np.array(data_df['class'])=='p')]

    
    
    # Plot
    plt.rcParams['font.family']='serif'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(2)
    
    plt.hist(periods,range=[0,20],bins=16,color='#0051a2',edgecolor='black')
    plt.xlabel('Period (d)')
    plt.ylabel('Frequency')
        
    plt.tick_params(axis='both',which='both',direction='in')
    plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    plt.minorticks_on()
    plt.subplots_adjust(wspace=0,hspace=0)
    plt.show()
    #plt.savefig(save_path,dpi=500,bbox_inches='tight')
    
sig_histogram('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/random/all.csv',save_path='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/random/periods_dist_c3.png')