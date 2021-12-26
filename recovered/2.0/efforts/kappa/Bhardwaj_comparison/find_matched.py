def find_matched(key,bhardwaj_file,results_file): 
    # Import(s)
    import numpy as np
    import pandas as pd
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from progress.bar import Bar
    
    # Action
    key_df = pd.read_csv(key)
    key_ids = np.array(key_df['ID'])
    key_ras = np.array(key_df['RA'])
    key_decs = np.array(key_df['DEC'])
    
    bhard_df = pd.read_csv(bhardwaj_file)
    bhard_ras = np.array(bhard_df['RA'])
    bhard_decs = np.array(bhard_df['DEC'])
    bhard_periods = np.array(bhard_df['PER'])
    
    shared_ids = np.array([])
    shared_b_ras = np.array([])
    shared_b_decs = np.array([])
    shared_periods = np.array([])
    
    bar = Bar('Processing...',max=len(key_ids))
    
    for key_id, key_ra, key_dec in zip(key_ids,key_ras,key_decs):
        c1 = SkyCoord(str(key_ra),str(key_dec),frame='fk5',unit=(u.deg,u.deg))
        for bhard_ra, bhard_dec, period in zip(bhard_ras,bhard_decs,bhard_periods):
            c2 = SkyCoord(str(bhard_ra),str(bhard_dec),frame='fk5',unit=(u.deg,u.deg))
            sep = c1.separation(c2).arcsecond
            if sep < 2:
                shared_ids = np.append(shared_ids,key_id)
                shared_b_ras = np.append(shared_b_ras,bhard_ra)
                shared_b_decs = np.append(shared_b_decs,bhard_dec)
                shared_periods = np.append(shared_periods,period)
                break
            
        bar.next()
            
    bar.finish()
 
    print(len(shared_ids))
    
    # Save shared ids to file
    zipped = list(zip(shared_ids,shared_b_ras,shared_b_decs,shared_periods))
    results_df = pd.DataFrame(data=zipped,columns=['ID','B_RA','B_DEC','B_PER'])
    results_df.to_csv(results_file,index=False)
            
    
        
    
key_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results_RA_AND_DEC.csv'
bhardwaj_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Bhardwaj_comparison/Bhardwaj_results.csv'
results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Bhardwaj_comparison/shared_results.csv'
find_matched(key=key_file,bhardwaj_file=bhardwaj_file,results_file=results_file)