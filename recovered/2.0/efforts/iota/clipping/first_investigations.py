def first_attempt(ids,data_temp,plot_path,sigma):
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from progress.bar import Bar

    # Action
    bar = Bar('Processing...',max=len(ids))
    for id in ids: 
        df = pd.read_csv(data_temp.replace('+++',id))
        mjd = np.array(df['mjd'])
        mag = np.array(df['mag'])
        
        mean_mag = np.mean(mag)
        std_mag = np.std(mag)

        # Plot
        plt.scatter(mjd,mag,color='#408ee0',s=2)
        plt.axhline(y=(mean_mag+sigma*std_mag),color='red',lw=2)
        plt.axhline(y=(mean_mag-sigma*std_mag),color='red',lw=2)
        plt.xlabel('Date (MJD)')
        plt.ylabel('Mag')
        plt.gca().invert_yaxis()
        plt.savefig(plot_path.replace('+++',id))
        plt.clf()
        plt.close()

        bar.next()

    bar.finish()

data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++_r.csv'
plot_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/clipping/5_sigma_test_plots/+++_r.png'
import pandas as pd
ids_list = list(pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/refined_good_ids.csv')['ID'])
first_attempt(ids=ids_list,data_temp=data_temp,plot_path=plot_temp,sigma=5)