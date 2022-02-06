def first():
    import pandas as pd
    import numpy as np
    import re 
    
    df1 = pd.read_csv('./referee_changes/new_names/old_new_names.csv')

    old_ids = list(df1['old_name'])
    new_ids = list(df1['preferred_name'])

    changed = 0
    for i in range(0,len(old_ids)):
        old = re.sub(' +', '_', old_ids[i])
        new = re.sub(' +', '_', new_ids[i])
        if old!=new: 
            changed+=1

    print('From LAH Changed:', changed)

    df2 = pd.read_csv('./results/meta_data.csv')
    old_ids = list(df2['ID'])
    new_ids = list(df2['preferred_name'])

    changed = 0
    for i, j in zip(old_ids, new_ids): 
        if i!=j:
            changed+=1

    print('Incorporated Changes:', changed)

    # both were 184

def second():
    import pandas as pd
    import numpy as np
    import re 
    
    df1 = pd.read_csv('./results/meta_data.csv')

    ids = np.array(df1['ID'])
    meta_preferred = np.array(df1['preferred_name'])

    df2 = pd.read_csv('./referee_changes/JANUARY_FINAL_RESULTS/JANUARY_FINAL_RESULTS.csv')
    ids_jan = np.array(df2['ID'])

    df3 = pd.read_csv('./referee_changes/new_names/old_new_names.csv')

    lah_preferred_ids = np.array(df3['preferred_name'])

    lah_preferred_ids = np.array([re.sub(' +', '_', i) for i in lah_preferred_ids])

    counter = 0
    for preferred in lah_preferred_ids: 
        if preferred not in meta_preferred: 
            print(preferred)
            counter+=1

    print(counter)

    # counter was 0


def third():
    import os
    import pandas as pd
    import numpy as np
    import re

    df3 = pd.read_csv('./referee_changes/new_names/old_new_names.csv')

    lah_preferred_ids = np.array(df3['preferred_name'])

    lah_preferred_ids = np.array([re.sub(' +', '_', i) for i in lah_preferred_ids])
    
    missing_counter = 0

    plot_ids = []
    for file in os.listdir(r"C:\Users\Research\Documents\GitHub\NAPYSOs\referee_changes\final_supp_plots\all_plots"):
        id = file.split('.')[0]
        plot_ids.append(id)
        
    for preferred in lah_preferred_ids: 
        if id not in plot_ids: 
            print(id)
            missing_counter+=1

    print(missing_counter)
    # none were missing! 


def fourth(): ## please work :praying_hands: 
    import pandas as pd
    import re
    
    df_meta = pd.read_csv('./results/meta_data.csv')
    meta_old = list(df_meta['ID'])
    meta_new = list(df_meta['preferred_name'])

    df_lah = pd.read_csv('./referee_changes/new_names/old_new_names.csv')
    lah_preferred = list(df_lah['preferred_name'])
    lah_old = [re.sub(' +', '_', i) for i in df_lah['old_name']]

    issue_counter = 0
    for lah_pref, lah_older in zip(lah_preferred, lah_old): 
        meta_index = meta_old.index(lah_older)
        if meta_new[meta_index]!=re.sub(' +', '_', lah_pref): 
            issue_counter+=1
            print(lah_pref, lah_older, meta_new[meta_index], meta_old[meta_index])

    print(issue_counter)
    # issues was 0! 

#first() --> issues: 0 (check how many lah changed against how many different in final results)
#second() --> issues: 0 (check how many lah changes were incorporated into final results)
#third() --> issues: 0 (check plots against lah new names)
#fourth() --> issues: 0 (check by matching indices based on value)