def copy_them(results_file, old_dir, new_dir):
    # Action
    import pandas as pd
    import numpy as np
    import os
    
    # Action
    slopes_df = pd.read_csv('/home/thaddaeus/June_6th.csv')
    df = pd.read_csv(results_file)
    df = df.join(slopes_df, on='ID')
    
    angle_errs = df['angle error']
    angles_mask = np.intersect1d(np.where(angle_errs<10), np.where(angle_errs>0))
    
    df = df.iloc[angles_mask]
    
    classes = np.array(df['primary_class'])
    classes_mask = np.where(classes=='p')
    
    df = df.iloc[classes_mask]
    
    print(df)
    
results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Misc./finalized_results/merged_results.csv'
old_dir = ''
new_dir = ''
copy_them(results_file, old_dir, new_dir)