def copy_plots(results_file,diagnostic_dir,one_temp,two_temp,three_temp):
    # Import(s)
    import numpy as np
    import pandas as pd
    import shutil
    import os
    from progress.bar import Bar
    
    # Action
    key_df = pd.read_csv(results_file)
    key_ids = np.array(key_df['filename'])
    key_classes = np.array(key_df['class'])
    
    one_indices = np.where(key_classes==1)
    one_ids = key_ids[one_indices]
    
    two_indices = np.where(key_classes==2)
    two_ids = key_ids[two_indices]
    
    three_indices = np.where(key_classes==3)
    three_ids = key_ids[three_indices]
    
    q_lm_filenames = os.listdir(diagnostic_dir)
    
    bar = Bar('Processing...',max=len(q_lm_filenames))
    
    for i in q_lm_filenames:
        i_list = i.split(':')
        name = i
        #name = i_list[2].replace('.png','')
        for one_id in one_ids: 
            if name==one_id: 
                shutil.copyfile(os.path.join(diagnostic_dir,i),os.path.join(one_temp,i))
                break
        
        for two_id in two_ids: 
            if name==two_id:
                shutil.copyfile(os.path.join(diagnostic_dir,i),os.path.join(two_temp,i))
                break
            
        for three_id in three_ids:
            if name == three_id: 
                shutil.copyfile(os.path.join(diagnostic_dir,i),os.path.join(three_temp,i))
                break
            
        bar.next()
    
    bar.finish()
    



results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Determing_periodicity/April_23/classes.csv'
diagnostic_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/plots/labeled_by_q'
good_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Determing_periodicity/April_23/plots/Periodic'
quest_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Determing_periodicity/April_23/plots/Questionable_Rejects'
reject_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Determing_periodicity/April_23/plots/Not_Periodic'

copy_plots(results_file,diagnostic_dir,good_temp,quest_temp,reject_temp)