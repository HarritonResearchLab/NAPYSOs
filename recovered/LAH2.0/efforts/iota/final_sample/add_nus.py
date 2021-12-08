def add_nus(key,data_temp):
    # Import(s)
    import numpy as np
    import pandas as pd

    # Definitions
    def sokolovskyNu(x,xerr):
        '''
        Best working version for calculating nu variabilty metric 
        from Sokolovsky et al. 2017. 
        '''
        
        #Import(s)
        import numpy as np
        
        #Action
        # [ (m-e)max - (m+e)min ] / [ (m-e)max + (m+e)min ]

        x = np.array(x)
        xerr = np.array(xerr)
        
        x_minus_xerr = x-xerr
        x_plus_xerr = x+xerr
        
        min_val = np.amin(x_plus_xerr)
        max_val = np.amax(x_minus_xerr)
        
        nu = (max_val-min_val)/(max_val+min_val)
        return nu

    # Action
    key_lines = []
    with open(key,'r') as f:
        for line in f:
            line = line.replace('\n','')
            key_lines.append(line)
    
    mod_key_lines = []
    for line_index, line in enumerate(key_lines): 
        if line_index == 0: 
            mod_key_lines.append(line)
        
        else: 
            line_list = line.split(',')
            id = line_list[0]
            data_file = data_temp.replace('+++',id)
            data_df = pd.read_csv(data_file)
            mags = np.array(data_df['mag'])
            magerrs = np.array(data_df['magerr'])
            nu = sokolovskyNu(x=mags,xerr=magerrs)
            line = line+','+str(nu)
            mod_key_lines.append(line)

    with open(key,'w') as f2:
        print('Erased') 

    with open(key,'a') as f3:    
        for line in mod_key_lines: 
            f3.write(line+'\n')

key = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results.csv'
data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
add_nus(key=key,data_temp=data_temp)