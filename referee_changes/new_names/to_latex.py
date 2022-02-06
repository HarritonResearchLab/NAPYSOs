import numpy as np
import pandas as pd


def sig_round(val, precision):
    str_per = str(round(val, 6))
    split_per = str_per.split('.')
    first_half = split_per[0]
    second_half = split_per[1]

    sig_per = first_half+'.'+second_half[0:precision-len(first_half)]

    return str(float(sig_per))

df = pd.read_csv('./results/meta_data.csv')
df = df.sort_values(by='preferred_name', ascending=False)

ids = np.array([i.replace('_', ' ') for i in df['preferred_name']])
ras = [str(round(i, 5)) for i in df['RA']]
decs = [str(round(i, 5)) for i in df['DEC']]
periods = list(df['period'].fillna('-'))

for index, i in enumerate(periods): 
    if i!='-': 
        periods[index] = sig_round(i,4)

qs = [str(round(i, 2)) for i in df['Q']]
ms = [str(round(i, 2)) for i in df['M']]
nus = [str(round(i, 4)) for i in df['nu']]
rs = [str(round(i, 2)) for i in df['r_mag']]
primaries = [i.upper() for i in df['primary_class']]
secondaries = [i.upper() for i in df['secondary_class'].fillna('-')]

angles = np.array(df['cmd_angle']).astype(object)
angle_errs = np.array(df['angle_error'])

for index, i in enumerate(angle_errs): 
    if angles[index] < 0: 
        angles[index] = str(round(angles[index]+180, 1))
    
    else: 
        angles[index]= str(round(angles[index], 1))

    if i<0 or i>10: 
        angles[index] = '-'

zipped = list(zip(ids, ras, decs, periods, 
                  qs, ms, nus, rs, primaries, secondaries, angles))

tex_path = './referee_changes/new_names_redux/table_2_body.tex'

latex_df = pd.DataFrame(zipped).to_latex(buf=tex_path, index=False)

lines = []
midrule_idx = 0
bottom_index = 0
line_index = 0

with open(tex_path, 'r') as f:
    for line in f:

        if '\midrule' in line: 
            midrule_idx = line_index+1
        
        if '\end{tabular}' in line: 
            bottom_index = line_index-1
            break
        
        lines.append(line)
        line_index+=1
        
lines = lines[midrule_idx:bottom_index]


with open(tex_path, 'w') as f: 
    for i in lines: 
        f.write(i)
