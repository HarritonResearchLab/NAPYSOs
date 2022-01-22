import pandas as pd
import numpy as np

df = pd.read_csv('referee_changes/JANUARY_FINAL_RESULTS/JANUARY_FINAL_RESULTS.csv')

primaries = np.array(df['primary_class'])
secondaries = np.array(df['secondary_class'])

classes_dict = {'p':'Periodic', 'mp':'Multi-periodic', 
                'qps':'Quasi-periodic symmetric', 
                'qpd':'Quasi-periodic dipper',
                'apd':'Aperiodic dipper', 'b':'Burster',
                'l':'Long timescale', 's':'Stochastic',
                'u':'Unclassifiable'}

len_one = len(primaries)
len_two = len(secondaries)

primary_percents = []
secondary_percents = []

for i, label in zip(classes_dict.keys(), classes_dict.values()):
    primary_percentage = round(len(np.where(primaries==i)[0])/len_one*100, 1)
    primary_percents.append(len(np.where(primaries==i)[0])/len_one*100)
    primary_percentage = str(primary_percentage)
    
    secondary_percentage = round(len(np.where(secondaries==i)[0])/len_two*100, 1)
    secondary_percents.append(len(np.where(secondaries==i)[0])/len_two*100)
    secondary_percentage = str(secondary_percentage)

    if secondary_percentage == '0.0': 
        secondary_percentage = '-'

    print(label + ' & ' + primary_percentage + ' & ' + secondary_percentage + r' \\')


print('\n****\n')

print('Total & ' + str(round(sum(primary_percents), 0)) + ' & '+str(round(sum(secondary_percents), 1))+r' \\')
