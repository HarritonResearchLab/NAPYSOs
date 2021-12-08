def remove_periodics(qm_file,periods_file,old_plot_dir,new_plot_dir,out_file):
    # Import(s)
    import numpy as np
    import pandas as pd
    import os
    import shutil 
    
    # Action
    
    qm_df = pd.read_csv(qm_file)
    qm_ids = np.array(qm_df['ID'])
    qm_ms = np.array(qm_df["M"])
    qm_qs = np.array(qm_df['Q'])
    qm_nus = np.array(qm_df['NU'])
    
    pers_df = pd.read_csv(periods_file)
    pers_ids = np.array(pers_df['filename'])
    pers_class = np.array(pers_df['class'])
    
    q_labeled_file_names = np.array(os.listdir(old_plot_dir))
    
    with open(out_file,'a') as f:
        f.write('ID,Q,M,NU,class'+'\n')
        for qm_id, m, q, nu in zip(qm_ids,qm_ms,qm_qs,qm_nus):
            for file_name, per_class in zip(pers_ids,pers_class):
                match = False
                if per_class == 1: 
                    if qm_id in file_name: 
                        match = True
                        f.write(str(qm_id)+','+str(q)+','+str(m)+','+str(nu)+','+'p'+'\n')
                        break
                        
            if match == False: 
                for q_file_name in q_labeled_file_names:
                    if qm_id in q_file_name: 
                        #shutil.copyfile(os.path.join(old_plot_dir,q_file_name),os.path.join(new_plot_dir,q_file_name))
                        break
                
                
    

qm_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results.csv'
periods_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Determing_periodicity/April_23/classes.csv'
old_plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/plots/labeled_by_q' 
new_plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/qm_classes_alpha/non_periodic_plots'
out_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/qm_classes_alpha/qm_classes.csv'
remove_periodics(qm_file,periods_file,old_plot_dir,new_plot_dir,out_file)
