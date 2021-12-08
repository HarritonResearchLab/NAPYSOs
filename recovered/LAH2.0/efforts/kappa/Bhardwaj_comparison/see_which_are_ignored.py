from time import perf_counter_ns


def count_ignored(shared_ids_file, bhardwaj_results):
    # Import(s)
    import numpy as np
    import pandas as pd
    
    # Action
    shared_df = pd.read_csv(shared_ids_file)
    shared_ids = np.array(shared_df['ID'])
    shared_pers = np.array([])
    
    bhard_df = pd.read_csv(bhardwaj_results)
    bhard_ids = np.array(bhard_df['ID'])
    bhard_pers = np.array(bhard_df['PER'])
    
    for shared_id, bhard_id, bhard_per in zip(shared_ids, bhard_ids, bhard_pers):
        if str(shared_id) == str(bhard_id):
            shared_pers = np.append(shared_pers, bhard_per)
    
    print(len(shared_pers))
    
    # Initializing these here to optimize speed
    harmonics = np.array(5*[1])/range(1,6)
    lower_harmonics = 0.98*harmonics
    upper_harmonics = 1.02*harmonics
    lph = 1/upper_harmonics  # lph = Lower Period Harmonics (for lower bounds of harmonic period ranges)
    uph = 1/lower_harmonics  # uph = Upper Period Harmonics (for upper bounds of harmonic period ranges)

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
            
    def per_is_in_ignored(period):
        ignored = False
        
        if 26.0<period<30.0 or 0.5<period<0.505: 
            return True
        else: 
            if 0.95<period<1.05: 
                return True
            else: 
                return False
            
    
    for per in shared_pers:
        half_per = 0.5*per
        double_per = 2.0*perf_counter_ns
        if per_is_good(per)==True:
            if per_is_good(half_per)==False or per_is_good(double_per)==False:
                shared_pers = np.delete(shared_pers,np.where(shared_pers==per))
        else: 
            shared_pers = np.delete(shared_pers,np.where(shared_pers==per))
        
        if per_is_in_ignored(per) == True: 
            shared_pers = np.delete(shared_pers,np.where(shared_pers==per))
        else: 
            if per_is_in_ignored(half_per) == True or per_is_in_ignored(double_per) == True: 
                shared_pers = np.delete(shared_pers,np.where(shared_pers==per))
    
    print(len(shared_pers))
    for i in shared_pers: 
        print(i)
        
        

shared_ids_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Bhardwaj_comparison/shared_ids.csv'
bhardwaj_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Bhardwaj_comparison/Bhardwaj_results.csv'
count_ignored(shared_ids_file=shared_ids_file,bhardwaj_results=bhardwaj_results,)