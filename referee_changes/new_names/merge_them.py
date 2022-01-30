import pandas as pd
import numpy as np
import re

current_df = pd.read_csv('./referee_changes/JANUARY_FINAL_RESULTS/full_results.csv')

new_names_df = pd.read_csv('./referee_changes/new_names/old_new_names.csv')

all_current_ids = np.array(current_df['ID'])
new_current_ids = np.array([re.sub(' +', '_', i) for i in new_names_df['old_name']])
new_ids = np.array([re.sub(' +', '_', i) for i in new_names_df['preferred_name']])


mixed_ids = []

print(len(new_ids))

for old_id in all_current_ids: 
    not_found = True
    for id_index, current_id in enumerate(new_current_ids): 
        if old_id ==current_id: 
            mixed_ids.append(new_ids[id_index])
            not_found = False
            break
    
    if not_found: 
        mixed_ids.append(old_id)

current_df['preferred_name']=np.array(mixed_ids)
current_df.to_csv('renamed.csv', index=False)

for old, true_old, final in zip(all_current_ids, current_df['ID'], mixed_ids):
    print(old,true_old, final)