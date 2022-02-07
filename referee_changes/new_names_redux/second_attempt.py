import pandas as pd
import numpy as np
import re

def first():
    lines = []
    with open('./referee_changes/new_names_redux/names_FHK_lahedits (1).csv', 'r') as f: 
        for line in f: 
            if line!='\n':
                lines.append(line)

    with open('./referee_changes/new_names_redux/names_FHK_lahedits (1).csv', 'w') as f: 
        print('okay')

    with open('./referee_changes/new_names_redux/names_FHK_lahedits (1).csv' , 'a') as f:     
        for index, line in enumerate(lines):
            f.write(line) 

def second(): 
    lah_df = pd.read_csv('./referee_changes/new_names_redux/names_FHK_lahedits (1).csv')

    lah_df['old_name'] = [re.sub(' +', '_', i) for i in lah_df['old_name']]
    lah_df['preferred_name'] = [re.sub(' +', '_', i) for i in lah_df['preferred_name']]
    
    meta_df = pd.read_csv('./results/meta_data.csv')


    meta_df = meta_df.drop(columns=['preferred_name'])
    lah_df = lah_df.drop(columns=['ra','dec','Simbad','TMASS','GaiaDR2','NameFHK'])


    final_meta = meta_df.merge(lah_df, left_on='ID', right_on='old_name')

    #final_meta.to_csv('for_lah.csv', index=False)

    final_df = final_meta.drop(columns=['old_name'])
    #final_df = final_df.rename(columns={'preferred_name':'ID'})
    
    final_df_cols = list(final_df)[:-1]
    final_df_cols.insert(0, 'preferred_name')

    final_df = final_df[final_df_cols]
    
    final_df.to_csv('./results/results.csv', index=False)

second()