def probe_df(results_file): 
    # Import(s)
    import numpy as np
    import pandas as pd
    
    # Action
    df = pd.read_csv(results_file)
    
    ids = np.array(df['ID'])
    qs = np.array(df['Q'])
    ms = np.array(df['M'])
    pers = np.array(df['PER1'])
    nus = np.array(df['NU'])
    
    for id, q, m, per, nu in zip(ids, qs, ms, pers, nus):
        if id == '2MASS_J20470481+4349114':
            print(q)
            print(m)
            print(per)
            print(nu)

results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/Modifying_LaTeX_Code/stuff_for_results_table/all.csv'
probe_df(results_file)