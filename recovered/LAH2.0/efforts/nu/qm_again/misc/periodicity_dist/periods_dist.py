def plot_rotation_periods(ref_file, plot_path):
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.ticker import AutoMinorLocator

    # Action
    df = pd.read_csv(ref_file)
    df = df.dropna()
    all_periods = np.array(df['PER'])
    primary_classes = np.array(df['primary_class'])
    good_periods = all_periods[np.where(primary_classes=='p')]
    #good_periods = all_periods[np.where(all_periods!=np.nan)]
    #good_periods = good_periods[np.logical_or(good_periods>1.05, good_periods<0.95)]
    
    #good_periods = np.sort(good_periods[np.where(good_periods>0.55)])
    
    # Plot
    plt.rcParams['figure.figsize']=(4,4)
    plt.rcParams['font.family'] = 'serif'
    
    # For all periods: 
    #plt.hist(good_periods,bins=20,color='#408ee0',edgecolor='black',range=(0.5,25))
    
    # For p only: 
    plt.hist(good_periods,color='#408ee0',edgecolor='black')
    
    plt.tick_params(axis='both',which='both',direction='in')
    plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
    plt.xlabel('Period (d)')
    plt.ylabel('Frequency')
    #plt.show()
    plt.savefig(plot_path, dpi=400, bbox_inches='tight')

# ACTUALLY RUN IT
ref_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Misc./finalized_results/merged_results.csv'
plot_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/plots for lah edits/periods_dist.png'
plot_rotation_periods(ref_file, plot_path)