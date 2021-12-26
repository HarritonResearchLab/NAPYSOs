import pandas as pd
import numpy as np

df = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv')
classes = list(set(df['primary_class']))

primary_classes = np.array(df['primary_class'])
secondary_classes = np.array(df['secondary_class'])

for qm_class in classes:
    num_class = len(primary_classes[np.where(primary_classes==qm_class)])
    print('PRIMARY: ', qm_class, ' num: ', num_class,  ' % ', 100*num_class/len(primary_classes))
    num_class = len(secondary_classes[np.where(secondary_classes==qm_class)])
    print('SECONDARY: ', qm_class, ' num: ', num_class,  ' % ', 100*num_class/len(secondary_classes))