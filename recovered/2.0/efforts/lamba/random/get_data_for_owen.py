def first_step(class_file):
    # Import(s)
    import pandas as pd
    import numpy as np
    
    # Action
    df = pd.read_csv(class_file)
    filenames = np.array(df['filename'])
    ids = np.array([])
    for filename in filenames: 
        ID = filename.split(':')[2].replace('.png','')
        ids = np.append(ids,ID)
    
    classes = np.array(df['class'])
    
    for ID, qm_class in zip(ids,classes):
        if (qm_class!='p') and (qm_class!='qps') and (qm_class!='qpd') and (qm_class!='apd') and (qm_class!='b') and (qm_class!='s') and (qm_class!='l') and (qm_class!='u'): 
            print(ID,qm_class)
            
            
#first_step('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/qm_classes_alpha/non_periodic_qm_classes.csv')

def second_step(results_file,data_temp,out_file):
    # Import(s)
    import numpy as np
    import pandas as pd
    
    # Action
    df = pd.read_csv(results_file)
    ids = np.array(df['ID'])
    mean_mags = np.array([])
    for ID in ids: 
        data_df = pd.read_csv(data_temp.replace('+++',ID))
        mags = np.array(data_df['mag'])
        mean_mag = np.mean(mags)
        mean_mags = np.append(mean_mags,mean_mag)
    
    zipped = list(zip(ids,mean_mags))
    out_df = pd.DataFrame(zipped,columns=['ID','mean_r'])
    out_df.to_csv(out_file,index=False)
    
#second_step('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/random/stuff_for_results_table/all_raw_results_RA_AND_DEC.csv',
#'/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv',
#'/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/random/stuff_for_results_table/mean_r_mags.csv')

def third_and_final(data_file,out_file):
    # Import(s)
    import numpy as np
    import pandas as pd
    
    # Action
    df = pd.read_csv(data_file)
    file_names = np.array(df['filename'])
    classes = np.array(df['class'])
    
    ids = np.array([])
    for filename in file_names: 
        ID = filename.split(':')[2].replace('.png','')
        ids = np.append(ids,ID)
        
    zipped = list(zip(ids,classes)) 
    out_df = pd.DataFrame(zipped,columns=['ID','qm_class'])
    out_df.to_csv(out_file,index=False)   
    
third_and_final('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/qm_classes_alpha/non_periodic_qm_classes.csv',
'/home/thaddaeus/FMU/HRL/LAH2.0/efforts/lamba/random/stuff_for_results_table/qm_nonperiodic_classes.csv')