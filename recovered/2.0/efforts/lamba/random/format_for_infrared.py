def format_data_files(input_file,output_file):
    # Import(s)
    import numpy as np
    import pandas as pd
    
    # Action 
    data_df = pd.read_csv(input_file)
    ids = np.array(data_df['ID'])
    ras = np.array(data_df['RA'])
    decs = np.array(data_df['DEC'])
    
    zipped = list(zip(ras,decs,ids))
    zipped_df = pd.DataFrame(zipped,columns=['RA','DEC','ID'])    
    zipped_df.to_csv(output_file,index=False)    
    
input_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/key.csv'
output_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/random/key_for_infrared.csv'
format_data_files(input_file,output_file)