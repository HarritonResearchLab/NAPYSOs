def examine_periodicity(results_file,confidence):
    #Import(s)
    import pandas as pd
    
    #Action
    df = pd.read_csv(results_file)
    ids = list(df['Name'])
    periods = list(df['Best Period (days)'])
    powers = list(df['Power'])
    faps = list(df[(str(confidence)+'% FAP')])

    zipped = list(zip(ids,periods,powers,faps))

    sig_periods = [] #periods that are significant 
    sig_ids = [] #ids of objects with significant periods
    for id, period, power, fap in zipped:
        if power > 2*fap:
            sig_ids.append(id)
            sig_periods.append(period)
    
    print(100*(len(sig_periods)/len(periods)))

    def plot_periods():
        #Import(s)
        import matplotlib.pyplot as plt

        #Action

        plt.hist(sig_periods,color='cornflowerblue',bins=80)
        #plt.xscale('log')
        plt.xlabel('Period')
        plt.ylabel('Frequency')
        plt.show()

    plot_periods()
def q_metric(data):
    #Import(s)
    from scipy import stats 
    import numpy as np
    
    #Action
    data = np.array(data)
    q_metric = (np.mean([stats.mstats.mquantiles(data,prob=0.9),stats.mstats.mquantiles(data,prob=0.1)])-np.median(data))/np.sqrt(((data-data.mean())**2).sum()/len(data))
    return q_metric

#examine_periodicity(results_file='ls_r_results.csv',confidence=99)  

def demo_plot():
    import matplotlib.pyplot as plt

    fig,axs = plt.subplots(2,3)


    axs[0,0].plot([],[],label='med g < 20.5: Slope < 2 vs Flux asym')
    axs[0,0].legend(loc='center',fontsize='xx-small')

    axs[0,1].plot([],[],label='med g < 20.5: Slope > 2 vs Flux asym')
    axs[0,1].legend(loc='center',fontsize='xx-small')

    axs[0,2].plot([],[],label='med g < 20.5: All Slopes vs Flux asym')
    axs[0,2].legend(loc='center',fontsize='xx-small')

    axs[1,0].plot([],[],label='All Slopes vs Flux asym')
    axs[1,0].legend(loc='center',fontsize='xx-small')

    axs[1,1].plot([],[],label='Slope > 2 vs Flux asym')
    axs[1,1].legend(loc='center',fontsize='xx-small')

    axs[1,2].plot([],[],label='All Slopes vs Flux asym')
    axs[1,2].legend(loc='center',fontsize='xx-small')

    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    
    plt.show()

demo_plot()

    