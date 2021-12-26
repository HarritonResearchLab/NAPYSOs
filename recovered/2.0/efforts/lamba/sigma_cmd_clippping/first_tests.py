def sigma_vs_percentile_plots(key,data_temp,plot_dir):
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    from progress.bar import Bar
    
    # Action
    key_df = pd.read_csv(key)
    ids = np.array(key_df['ID'])
    
    bar = Bar('Processing...',max=len(ids))
    
    plt.style.use('seaborn-darkgrid')
    plt.rcParams['font.family']='serif'
    
    for id in ids: 
        data = pd.read_csv(data_temp.replace('+++', id))
        g_r = np.array(data['g-r_val'])
        g = np.array(data['g_mag'])
        
        mean_color = np.mean(g_r)
        std_color = np.std(g_r)
        
        mean_g = np.mean(g)
        std_g = np.std(g)
        
        left_edge = np.percentile(g_r,2.5)
        right_edge = np.percentile(g_r,97.5)
        
        top_edge = np.percentile(g,97.5)
        bottom_edge = np.percentile(g,2.5)
        
        plt.scatter(g_r,g,color='#408ee0',edgecolor='black')
        
        # Percentile Edges
        plt.axvline(x=left_edge,color='red', label='Middle 95%')
        plt.axvline(x=right_edge,color='red')
        plt.axhline(y=top_edge,color='red')
        plt.axhline(y=bottom_edge,color='red')
        
        # 5 Sigma Edges
        plt.axvline(x=mean_color-4*std_color,color='blue',label='4 '+r'$\sigma$')
        plt.axvline(x=mean_color+4*std_color,color='blue')
        #plt.axhline(y=mean_g-5*std_g,color='blue')
        #plt.axhline(y=mean_g+5*std_g,color='blue')
        
        plt.legend()
        
        plot_path = os.path.join(plot_dir,(id+'.png'))
        plt.savefig(plot_path)
        plt.clf()
        plt.close()
        
        bar.next()
    
    bar.finish()
        

key = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/key.csv'
data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/cmd/+++.csv'
plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/sigma_cmd_clippping/plots'
sigma_vs_percentile_plots(key,data_temp,plot_dir)