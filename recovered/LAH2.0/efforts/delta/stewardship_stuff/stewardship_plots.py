def mag_vs_var(key,path_temp,band):
    # Import(s)
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    from personalastropy.ysospy.variability import sokolovskyNu
    from progress.bar import ShadyBar

    # Action

    key_df = pd.read_csv(key)
    ids = list(key_df['ID'])

    nus = [] # Sok. nu metrics
    stds = [] # stds. of mag dists. 
    mean_mags = [] # mean mags of mag dists. 

    bar = ShadyBar('Processing...',max=len(ids))

    for id in ids:
        data_file = path_temp.replace('+++',id)
        data_df = pd.read_csv(data_file)
        
        mags = np.array(data_df['mag'])

        if len(mags) > 10: 
            magerrs = np.array(data_df['magerr'])

            std = np.std(mags,ddof=1)
            stds.append(float(std))

            nu = sokolovskyNu(mags,magerrs)
            nus.append(float(nu))

            mean_mag = np.mean(mags)
            mean_mags.append(float(mean_mag))
        
        bar.next()

    bar.finish()

    # Plot
    plt.rcParams['font.family'] = 'serif'
    fig, axs = plt.subplots(2,1)
    
    axs[0].scatter(mean_mags,nus,color='cornflowerblue',marker='.',s=2)
    axs[0].set_ylabel('Sokolovsky '+r'$\nu$')

    axs[1].scatter(mean_mags,stds,color='cornflowerblue',marker='.',s=2)
    axs[1].set_ylabel(r'$\sigma$')
    axs[1].set_xlabel('Mean Mag')

    fig.suptitle(band+' Variability',fontsize='small')

    #plt.subplots_adjust(hspace=0)

    plt.show()

def q_vs_p(ls_results,qm_results):
    # Import(s)
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    # Action

    ls_df = pd.read_csv(ls_results)
    ls_ids = list(ls_df['ID'])
    ls_pers = list(ls_df['PER'])
    ls_faps = list(ls_df['99FAP'])
    ls_pows = list(ls_df['POW'])

    sig_ls_ids = []
    sig_ls_pers = []
    sig_ls_pows = []
    
    for i, id in enumerate(ls_ids):
        pow = ls_pows[i]
        fap = ls_faps[i]
        per = ls_pers[i]
        if pow>fap and pow>0.1: 
            sig_ls_ids.append(id)
            sig_ls_pers.append(per)
            sig_ls_pows.append(pow)

    final_pows = []
    final_pers = []
    qs = []

    qm_df = pd.read_csv(qm_results)
    qm_ids = list(qm_df['ID'])
    initial_qs = list(qm_df['Q'])

    for i, id in enumerate(sig_ls_ids):
        for qm_id in qm_ids:
            if id==qm_id:
                qm_index = qm_ids.index(qm_id)
                q = initial_qs[qm_index]
                qs.append(float(q))
                final_pows.append(sig_ls_pows[i])
                final_pers.append(sig_ls_pers[i])
                print(id,str(q))
                
    
    final_pows = np.array(final_pows)
    final_pers = np.array(final_pers)
    qs = np.array(qs)

    print(len(qs[np.where(qs>0.15)]))

    # Plot
    plt.rcParams['font.family'] = 'serif'
    scaled_pows = (10*(final_pows**2))+1
    scaled_pows = np.around(scaled_pows,decimals=2)
    plt.scatter(final_pers,qs,s=scaled_pows)
    plt.xlabel('Period (d)')
    plt.ylabel('Quasi-Periodicity (Q)')
    plt.show()

qm_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/gamma/testing_qm/q_vs_m_results.csv'
ls_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/delta/stewardship_stuff/redone_ls_again/report_file.csv'
q_vs_p(ls_results=ls_file,qm_results=qm_file)
