
from matplotlib.pyplot import minorticks_on


def comparison_plot(results_file, data_temp):
    # Import(s)
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    
    # Definitions
    def bredall_m(mags):
        # Import(s)
        import numpy as np

        # Action
        (tenth, ninetieth) = np.percentile(mags, [10, 90])
        mean = np.mean(mags[np.logical_or(mags>ninetieth, mags<tenth)])

        m = (mean-np.median(mags))/np.std(mags)
        return m
    
    # Action
    df = pd.read_csv(results_file)
    ids = np.array(df['ID'])
    our_ms = np.array(df['M'])
    
    correct_ms = np.array([])
    
    for id in ids: 
        data_file = pd.read_csv(data_temp.replace('+++', id))
        correct_m = bredall_m(data_file['mag'])
        correct_ms = np.append(correct_ms, correct_m)
        
    # Make plot
    
    fig, axs = plt.subplots(1, 2, figsize=(8, 3))
    plt.style.use('ggplot')
    
    axs[0].scatter(correct_ms, our_ms)
    axs[0].set(xlabel='Corrected M', ylabel='Current M')
    axs[0].axline((-1, -1), slope=1, label='1:1', color='gray', ls='--')
    
    
    axs[1].scatter(correct_ms, correct_ms-our_ms)
    axs[1].set(xlabel='Corrected M', ylabel='Corrected M - Current M')
    axs[1].axline((-1, -1), slope=1, label='1:1', color='gray', ls='--')
    axs[1].legend()
    
    for ax in [axs[0], axs[1]]:
        ax.minorticks_on()
        ax.legend()
        ax.axhline(y=0.25)
        ax.axhline(y=-0.25)
        ax.axvline(x=0.25)
        ax.axvline(x=-0.25)
    
    plt.subplots_adjust(wspace=0.5)
    #plt.show()
    plt.savefig('m_comparison.png', bbox_inches='tight')
        
results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv'
data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/light_curves/+++_r.csv'
comparison_plot(results_file, data_temp)