import pandas as pd
import numpy as np
import re

current_df = pd.read_csv('./referee_changes/JANUARY_FINAL_RESULTS/full_results.csv')

new_names_df = pd.read_csv('./referee_changes/new_names/old_new_names.csv')

all_current_ids = np.array(current_df['ID'])
new_current_ids = np.array([re.sub(' +', '_', i) for i in new_names_df['old_name']])
new_ids = np.array([re.sub(' +', '_', i) for i in new_names_df['preferred_name']])

mixed_ids = []

for current_id in all_current_ids: 
    if current_id in new_current_ids: 
        index = np.where(new_current_ids==current_id)[0]
        new_id = new_ids[index][0]


        mixed_ids.append(new_ids[index][0])

    else: 
        mixed_ids.append(current_id)

current_df['preferred_name']=np.array(mixed_ids)
current_df.to_csv('renamed.csv', index=False)

