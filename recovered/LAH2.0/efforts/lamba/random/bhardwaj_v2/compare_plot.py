def comparison_plot(our_results,bhards_results,old_plots_dir,new_plots_dir):
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import os 
    import shutil 
    from matplotlib.ticker import AutoMinorLocator
    
    # Action
    our_results = pd.read_csv(our_results)
    bhards = pd.read_csv(bhards_results)
    
    bhard_ids = np.array(bhards['ID'])
    our_ids = np.array(our_results['ID'])
    
    bhard_periods = np.array(bhards['B_PER'])
    our_periods = np.array(our_results['PER1'])
    second_pers = np.array(our_results['PER2'])
    third_pers = np.array(our_results['PER3'])
    
    faps = np.array(our_results['99%_FAP'])
    first_powers = np.array(our_results['POW1'])
    second_powers = np.array(our_results['POW2'])
    
    
    shared_b_periods = np.array([])
    shared_b_second_periods = np.array([])
    shared_our_periods = np.array([])
    shared_second_periods = np.array([])
    
    periodogram_file_names = np.array(os.listdir(old_plots_dir))
    periodogram_object_ids = np.array([])
    for file_name in periodogram_file_names:
        file_list = file_name.split(':')
        if file_list[0] != '.DS_Store' and file_list[0]!='results.csv':
            object_name = file_list[2]
            object_name = object_name.replace('.png','')
            periodogram_object_ids = np.append(periodogram_object_ids,object_name)
    
    for b_id, b_per in zip(bhard_ids,bhard_periods):
        for our_id, our_per, second_per, first_pow, second_pow, fap in zip(our_ids,our_periods,second_pers,first_powers, second_powers, faps):
            if b_id == our_id: 
                if first_pow > fap: 
                    shared_b_periods = np.append(shared_b_periods,b_per)
                    shared_our_periods = np.append(shared_our_periods,our_per)
                    
                    div = b_per/our_per
                    if div>1.2 or div<0.8: 
                        for id, file_name in zip(periodogram_object_ids,periodogram_file_names):
                            if our_id == id: 
                                shutil.copyfile(os.path.join(old_plots_dir,file_name),os.path.join(new_plots_dir,file_name))
        
                if second_pow > fap: 
                    shared_b_second_periods = np.append(shared_b_second_periods,b_per)
                    shared_second_periods = np.append(shared_second_periods,second_per)
                
                break

    # Plot
    plt.rcParams['font.family']='serif'
    plt.gca().tick_params(axis='both',which='both',direction='in')
    
    
    
    one_to_one = np.linspace(min(shared_our_periods),15,75)
    plt.plot(one_to_one,one_to_one,color='black',label='1:1',ls='--')
    plt.scatter(shared_our_periods,shared_b_periods,color='#408ee0',edgecolor='black', label='Best')
    plt.scatter(shared_second_periods,shared_b_second_periods,color='indianred',edgecolor='black',s=12,label='Second Best')
    #plt.scatter(shared_third_periods,bhard_periods,color='yellow',edgecolor='black',s=10,label='Third Best')
    plt.xlabel('Our Period')
    plt.ylabel('Bhardwaj et al. 2019 Period')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper right',fontsize=7,fancybox=False,edgecolor='black',shadow=False)
    #plt.axes().yaxis.set_minor_locator(AutoMinorLocator())
    #plt.axes().xaxis.set_minor_locator(AutoMinorLocator())
    plt.tick_params(axis='both',which='both',direction='in')
    plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
    plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    plt.savefig('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/random/bhardwaj_v2/log_log_comparison.png',dpi=500)
    
bhards_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Bhardwaj_comparison/shared_results.csv'
our_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/period_dists/Periodogram Plots_r/Periodogram Plots_r/results.csv'
old_plots_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/period_dists/Periodogram Plots_r/Periodogram Plots_r'
new_plots_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/random/bhardwaj_v2/outside_one_to_one'

comparison_plot(our_results,bhards_results,old_plots_dir,new_plots_dir)