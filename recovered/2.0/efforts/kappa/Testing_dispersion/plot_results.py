def plot_q_vs_quot(results_file):
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    
    # Action 
    results_df = pd.read_csv(results_file)
    
    ids = np.array(results_df['ID'])
    qs = np.array(results_df['Q'])
    quots = np.array(results_df['QUOTIENT']) 
    pows = np.array(results_df['POW'])
    faps = np.array(results_df['FAP'])
    
    sig_qs = np.array([])
    sig_quots = np.array([])
    
    for id, q, quot, pow, fap in zip(ids, qs, quots,pows,faps):
        if pow>fap: 
            sig_qs = np.append(sig_qs,q)
            sig_quots = np.append(sig_quots,quot)
            if q<0.60 and 1.25>quot>0.7: 
                print(str(id)+','+str(q)+','+str(quot))
    
    # Plot
    plt.rcParams['font.family'] = 'serif'
    fig, axs = plt.subplots(2,2)
    
    axs[0,0].scatter(qs,quots,color='#408ee0',edgecolor='black')
    axs[0,0].set_xlim(left=0,right=1)
    #axs[0,0].set_ylim(bottom=-)#,top=5)
    axs[0,0].tick_params(axis='both',which='both',direction='in')
    axs[0,0].tick_params(which='both',bottom=True,top=True,left=True,right=True)
    axs[0,0].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
    axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
    axs[0,0].set_xlabel('Q')
    axs[0,0].set_ylabel('Variance Quotient')
    axs[0,0].set_title('Variable Objects',fontsize='medium')
    
    axs[0,1].scatter(sig_qs,sig_quots,color='#408ee0',edgecolor='black')
    axs[0,1].set_xlim(left=0,right=1)
    #axs[0,1].set_ylim(bottom=0)#,top=5)
    axs[0,1].tick_params(axis='both',which='both',direction='in')
    axs[0,1].tick_params(which='both',bottom=True,top=True,left=True,right=True)
    axs[0,1].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
    axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
    axs[0,1].set_xlabel('Q')
    axs[0,1].set_ylabel('Variance Quotient')
    axs[0,1].set_title('Power > 99% FAP Objects',fontsize='medium')
    
    axs[1,0].axis('off')
    axs[1,1].axis('off')
    
    plt.subplots_adjust(wspace=0.5)
    plt.show()    
    
plot_q_vs_quot(results_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Testing_dispersion/results.csv')