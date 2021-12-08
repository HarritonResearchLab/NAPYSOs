def proto_merg(results, odr_results):
    # Imports(s)
    import pandas as pd
    
    # Action
    df1 = pd.read_csv(results)
    df2 = pd.read_csv(odr_results)
    
    df3 = df1.merge(df2, on='ID')
    df3.to_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv', index=False)
    
def format_for_table_2(merged_results, out_file):
    # Import(s)
    import numpy as np
    import pandas as pd
    from sigfig import round as sig_round
    
    # Action
    df = pd.read_csv(merged_results)
    df_length = len(df.index)
    
    formatted_lines = []
    
    # id, ra, dec, per, q, m, nu, r_mean, primary_class, secondary class, 
    # slope angle, angle error
    
    ids = np.array([i.replace('_', ' ') for i in df['ID']])
    ras = np.array(df["RA"])
    decs = np.array(df["DEC"])
    rmeans = np.array(df["rMean"])
    p_classes = np.array([str(i).upper() for i in df["primary_class"]])
    s_classes = df['secondary_class'].fillna('-')
    s_classes = np.array([str(i).upper() for i in s_classes])
    nus = np.array(df["NU"])
    periods = np.array(df["PER"].fillna('-'))
    qs = np.array(df["Q"])
    ms = np.array(df["M"])
    slope_angles = np.array(df['angle'])
    angle_errors = np.array(df['angle error'])
    
    for id, ra, dec, per, q, m, nu, r_mean, p_class, s_class, slope_angle, angle_error in zip(ids,
    ras, decs, periods, qs, ms, nus, rmeans, p_classes, s_classes, slope_angles, angle_errors):
        ra = round(ra, 5)
        dec = round(dec, 5)
        if type(per) != str: 
            per = round(per, 3)
        q = round(q, 2)
        m = round(m, 2)
        nu = round(nu, 4)
        r_mean = sig_round(r_mean, 4)
        
        if 0 < angle_error < 10: 
            slope_angle = str(round(slope_angle, 1))
        else: 
            slope_angle = '-'
        
        out_str = id +' & ' + str(ra) + ' & ' + str(dec) + ' & ' + str(per)
        out_str += ' & ' + str(q) + ' & ' + str(m) + ' & ' + str(nu) + ' & '
        out_str += str(r_mean) + ' & ' + p_class + ' & ' + s_class + ' & '
        out_str += slope_angle + r' \\' + '\n'
        
        formatted_lines.append(out_str)
        
    with open(out_file, 'a') as f:
        for line in formatted_lines: 
            f.write(line)
            
merged_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv'
out_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/formatted.txt'

format_for_table_2(merged_results, out_file)
        
    
        
    
    
        