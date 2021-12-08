def enumerate_classes(results_file):
    # Import(s)
    import pandas as pd
    import numpy as np
    
    # Action
    df = pd.read_csv(results_file)
    primary = np.array(df['primary_class'])
    secondary = np.array(df['secondary_class'])
    
    qm_classes = ['p', 'mp', 'qps', 'qpd', 'apd', 'b', 'l', 's', 'u']
    
    length = len(primary)
    
    print('PERCENTAGES')
    for qm_class in qm_classes: 
        primary_percent = 100*round(np.where(primary==qm_class)[0].size/length, 3)
        secondary_percent = 100*round(np.where(secondary==qm_class)[0].size/length, 3)
        print(qm_class.upper(), primary_percent, secondary_percent)
        
    print('\nCounts')
    for qm_class in qm_classes: 
        print(qm_class.upper(), np.where(primary==qm_class)[0].size)
        
            
enumerate_classes('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/september21/SEPTEMBER_FINAL_RESULTS/merged_results.csv')
    