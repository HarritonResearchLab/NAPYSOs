def get_oid(data_temp,id):
    # Import(s)
    from astropy.io import ascii
    
    # Action
    data_table = ascii.read(data_temp.replace('+++',id),format='ipac',delimiter='|')
    return data_table#np.array(data_table['programid'])

# Import(s)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Action
key_df = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/key.csv')
ids = list(key_df['ID'])

prog_ids = np.array([])

for id in ids: 
    prog_ids = np.append(prog_ids,get_oid('/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++_r.tbl',id))
     
plt.style.use('seaborn-darkgrid')  
plt.rcParams['font.family']='serif'
plt.hist(prog_ids)
plt.xlabel("Program ID")
plt.ylabel("Frequency")
plt.yscale("log")
plt.title('Distribution of programIDs',fontsize='medium')
plt.xticks(ticks=[1,2,3])
plt.show()
