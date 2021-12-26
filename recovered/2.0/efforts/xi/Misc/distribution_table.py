def make_distribution_table(results_file):
    # Import(s)
    import numpy as np
    import pandas as pd
    
    # Execution
    results_df = pd.read_csv(results_file)    
    primary_classes = np.array(results_df['primary_class'])
    results_df = results_df.dropna()
    secondary_classes = np.array(results_df['secondary_class'])
    
    classes = np.array(['p', 'mp', 'qps', 'qpd', 'apd', 'b', 'l', 's', 'u'])
    
    primary_length = len(primary_classes)
    secondary_length = len(secondary_classes)
    
    primary_sum = 0
    secondary_sum = 0
    for qm_class in classes: 
        print(qm_class.upper())
        
        primary_percentage = round(100*len(primary_classes[np.where(primary_classes==qm_class)])/primary_length, 1)
        secondary_percentage  = round(100*len(secondary_classes[np.where(secondary_classes==qm_class)])/primary_length, 1)
        
        primary_sum = primary_sum + primary_percentage
        secondary_sum = secondary_sum + secondary_percentage
            
results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/routine_results/finalized_results/merged_results.csv'
make_distribution_table(results_file)