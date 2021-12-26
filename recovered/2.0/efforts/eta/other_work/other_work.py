from numpy.core.defchararray import upper


def create_mags_vs_var_plot(key,data_path):
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    from personalastropy.ysospy.variability import sokolovskyNu

    # Action
    ids = list(pd.read_csv(key)['ID'])
    mean_mags = []
    nus = []
    stds = []
    for id in ids: 
        data_file = pd.read_csv(data_path.replace('+++',id))
        mags = np.array(data_file['mag'])
        magerrs = np.array(data_file['magerr'])
        mean_mag = np.mean(mags)
        nu = sokolovskyNu(x=mags,xerr=magerrs)
        std = np.std(mags)
        mean_mags.append(mean_mag)
        nus.append(nu)
        stds.append(std)
    
    mean_mags = np.array(mean_mags)
    nus = np.array(nus)
    stds = np.array(stds)

    min_mean = int(np.floor(min(mean_mags)))
    max_mean = int(np.ceil(max(mean_mags)))

    lower_bounds = range(min_mean,max_mean)

    # Quick inizialization of plot
    plt.rcParams['font.family'] = 'serif'
    plt.scatter(mean_mags,stds,s=1.5)    
    plt.xlabel('Mean r Mag')
    plt.ylabel('r '+r'$\nu$')


    for base_lower_bound in lower_bounds: 
        lower_bounds = [base_lower_bound,base_lower_bound+0.5]
        for lower_bound in lower_bounds: 
            upper_bound = lower_bound + 0.5
            mask = np.logical_and(mean_mags<upper_bound,mean_mags>=lower_bound)
            masked_nus = nus[mask]
            masked_stds = stds[mask]
            nu_percentiles = list(np.percentile(masked_nus,q=[10,50,90]))
            std_percentiles = list(np.percentile(masked_stds,q=[10,50,90]))
            
            colors = ['blue','orange','purple']
            for p_index, percentile in enumerate(nu_percentiles):
                xs = np.linspace(lower_bound,upper_bound,3)
                ys = np.array([percentile,percentile,percentile])
                plt.plot(xs,ys,color=colors[p_index],lw=2,ls='--')    
    plt.show()

dp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++_r.csv'
key_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/key.csv'
create_mags_vs_var_plot(key=key_path,data_path=dp)
