def investigate(key,data_file): 
    # Import(s)
    import numpy as np
    import pandas as pd
    from astropy.io import ascii
    
    # Action

    key = pd.read_csv(key)
    ids = list(key['ID'])

    for id in ids: 
        data_file = data_file.replace('+++',id)
        r = ascii.read(data_file,format='ipac',delimiter='|')
        catflags = np.array(r['catflags'])
        print(catflags[np.where(catflags==32768)])

investigate('/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/key.csv','/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++_r.tbl')