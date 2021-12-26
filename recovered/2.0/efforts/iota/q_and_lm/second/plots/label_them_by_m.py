def copy_and_rename(key,old_dir,new_dir):
    # Import(s)
    import numpy as np
    import pandas as pd
    import os
    import shutil
    from progress.bar import Bar

    # Action

    key_df = pd.read_csv(key)
    ids = list(key_df['ID'])
    ms = list(key_df['M'])

    bar = Bar('Processing...',max=len(ids))

    for id, m in list(zip(ids,ms)): 
        for file_name in os.listdir(old_dir):
            if id in file_name: 
                m_str = str(m)
                m_list = m_str.split('.')
                decimals = m_list[1]
                decimals = decimals[0:3]
                m_str = m_list[0]+'.'+decimals
                new_file_path = os.path.join(new_dir,(m_str+':ID:'+id+'.png'))
                old_file_path = os.path.join(old_dir,file_name)
                shutil.copyfile(old_file_path,new_file_path)

                bar.next()
            
    bar.finish()
            

# RUN IT
key = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results.csv'
old_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/plots/labeled_by_q'
new_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/plots/labeled_by_m'

copy_and_rename(key=key,old_dir=old_dir,new_dir=new_dir)