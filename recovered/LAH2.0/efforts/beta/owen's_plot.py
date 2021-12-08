def period_vs_gr(key,csv_temp):
    #Import(s)
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    #Action
    ls_df = pd.read_csv(key)
    ids = list(ls_df['ID'])
    periods = list(ls_df['Period'])
    powers = list(ls_df['Power'])
    faps = list(ls_df['99.9% FAP'])

    

    good_periods = []
    mean_grs = []
    
    for id, period, power, fap in zip(ids,periods,powers,faps):
        if power > fap:
            good_periods.append(period)
            cmd_df = pd.read_csv(csv_temp.replace('+++',id))
            g_minus_rs = np.array(cmd_df['g-r'])
            mean_gr = np.mean(g_minus_rs)
            mean_grs.append(mean_gr)
    


    #Plot
    plt.rcParams['font.family'] = 'serif'
    fig, axs = plt.subplots(2,1)

    axs[0].scatter(good_periods,mean_grs,color='cornflowerblue',marker='.')
    axs[0].set_xlim(0,50)

    axs[1].scatter(good_periods,mean_grs,color='cornflowerblue',marker='.')

    for i in range(2):
        axs[i].set_xlabel('Period')
        axs[i].set_ylabel('Mean g-r')

    plt.subplots_adjust(hspace=0.4)
    plt.show()




period_vs_gr(key='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/ls_results.csv',csv_temp='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/cmd/+++.csv')