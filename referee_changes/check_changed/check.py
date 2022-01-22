import pandas as pd 
import numpy as np

df = pd.read_csv('./referee_changes/lah_vetting_results/ref_files/lah_changes_incorporated.csv')

changed_mask = np.array(df['per_changed'].fillna(False))
changed_ids = np.array(df['ID'])[changed_mask]

with open('referee_changes/check_changed/raw_tex.txt', 'r') as f:
    for line in f: 
        for id in changed_ids: 
            if id in line: 
                print(line)