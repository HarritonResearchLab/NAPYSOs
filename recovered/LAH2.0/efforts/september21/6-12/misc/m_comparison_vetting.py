def get_m(id, data_temp):
    # Import(s)
    import pandas as pd
    import numpy as np
    
    # Action
    def bredall_m(mags):
        # Import(s)
        import numpy as np

        # Action
        (tenth, ninetieth) = np.percentile(mags, [10, 90])
        mean = np.mean(mags[np.logical_or(mags>ninetieth, mags<tenth)])

        m = (mean-np.median(mags))/np.std(mags)
        return m
    
    data_file = pd.read_csv(data_temp.replace('+++', id))
    correct_m = bredall_m(data_file['mag'])

    print(id, correct_m)
    
def return_sorted(results_file):
    # Import(s)
    import pandas as pd
    import numpy as np
    
    # Action
    df = pd.read_csv(results_file)
    ms = np.array(df['M'])
    ids = np.array(df['ID'])
    sorted_ms = ms[-15:]
    max_m = np.max(ms)
    print(ids[np.where(ms==max_m)], max_m)
    
        
results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv'
data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/light_curves/+++_r.csv'

get_m('GDR1_2162947420052773760', data_temp)