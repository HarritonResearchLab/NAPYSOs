def locate_missing_ones(key,g_path,r_path): 
    # Import(s)
    import pandas as pd
    import os

    # Action
    df = pd.read_csv(key)

    ids = list(df['ID'])
    sum = 0
    for id in ids: 
        g_file = g_path.replace('+++',id)
        r_file = r_path.replace('+++',id)

        if os.path.exists(g_file)==False or os.path.exists(r_file)==False: 
            print(id)

        sum = sum+1
    print(sum)
data_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_*.csv'
locate_missing_ones(key='/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/key.csv',g_path=data_path.replace('*','g'),r_path=data_path.replace('*','r'))
