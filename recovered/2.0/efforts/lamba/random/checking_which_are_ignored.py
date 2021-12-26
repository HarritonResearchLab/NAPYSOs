def check_them(data_file):
    # Imports(s)
    import pandas as pd
    import numpy as np
    
    # Action
    data_df = pd.read_csv(data_file)
    periods = np.array(data_df['B_PER'])
    
    harmonics = np.array(5*[1])/range(1,6)
    lower_harmonics = 0.98*harmonics
    upper_harmonics = 1.02*harmonics
    lph = 1/upper_harmonics # lph = Lower Period Harmonics (for lower bounds of harmonic period ranges)
    uph = 1/lower_harmonics # uph = Upper Period Harmonics (for upper bounds of harmonic period ranges)

    def per_is_good(period):
        harmonic = False
        for harmonic_idx in range(0,len(harmonics)):
            if lph[harmonic_idx] < period < uph[harmonic_idx]:
                harmonic = True
                break
        if harmonic==True: 
            return False
        else: 
            return True
            
    ignored_periods = np.array([])
    
    for per in periods: 
        if per_is_good(per)==False or 26<per<30 or 0.5<per<0.505: 
            ignored_periods = np.append(ignored_periods,per)
            
    print(ignored_periods)
    print(len(ignored_periods))

data_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Bhardwaj_comparison/shared_results.csv'
check_them(data_file)