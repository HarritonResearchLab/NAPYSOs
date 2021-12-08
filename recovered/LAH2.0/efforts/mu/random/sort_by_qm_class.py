from matplotlib.pyplot import plot


def sort_plots(results_file,old_plots_dir,new_plot_dir):
    # Import(s)
    import os
    import shutil
    import numpy as np
    import pandas as pd 
    
    # Action
    
    results_df = pd.read_csv(results_file)
    classes = np.array(results_df['class'])
    ids = np.array(results_df['ID'])
    
    old_ids = np.array([])
    old_filenames = np.array(os.listdir(old_plots_dir))
    for i in old_filenames: 
        if i != ".DS_Store": 
            old_ids = np.append(old_ids, i.split(':')[2].replace('.png',''))
    
    
    for qm_class in ['apd','qpd','qps','l','b','s','p','u','mp']:
        class_dir = os.path.join(new_plot_dir,qm_class)
        os.mkdir(class_dir)
        class_mask = np.where(classes==qm_class)
        masked_ids = ids[class_mask]
        for masked_id in masked_ids:
            for old_id, filename in zip(old_ids, old_filenames): 
                if masked_id == old_id: 
                    old_path = os.path.join(old_plots_dir,filename)
                    new_path = os.path.join(class_dir,filename)
                    shutil.copyfile(old_path,new_path)
                    break                     
            
    

results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/random/all.csv'
old_plots_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/random/Periodogram Plots_r'
new_plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/random/sorted_by_qm_class_new_format'
sort_plots(results_file,old_plots_dir,new_plot_dir)