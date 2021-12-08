#### NEED TO MAKE IT SUCH THAT ONLY SIG PERIODS ARE REPORTED!
### IN future, make new column for sig_period (Y/N) for periodic? 
# IDK doesn't seem nessesary  

# ## change per1 column in the all.csv file to just the 
# significant period for the object (which may be per2)

def compile():
    # Import(s)
    import pandas as pd
    import numpy as np
    
    # Action
    df_1 = pd.read_csv("/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/Modifying_LaTeX_Code/stuff_for_results_table/all_raw_results_RA_AND_DEC.csv")
    df_2 = pd.read_csv("/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/Modifying_LaTeX_Code/stuff_for_results_table/mean_r_mags.csv")
    df_3 = pd.read_csv("/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/Modifying_LaTeX_Code/stuff_for_results_table/qm_nonperiodic_classes.csv")
    df_4 = pd.read_csv("/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/Modifying_LaTeX_Code/stuff_for_results_table/qm_periodic_classes.csv")

    df_1.set_index("ID", inplace= True)
    df_2.set_index("ID", inplace= True)
    df_3.set_index("ID", inplace= True)
    df_4.set_index("ID", inplace= True)

    classes = df_4["class"]
    classes = classes.append(df_3["qm_class"])
    classes= pd.DataFrame(classes, columns = ["class"])

    total = df_1.join(df_2,how = "outer" )
    total = total.join(classes, how = "outer")
    total.to_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/Modifying_LaTeX_Code/stuff_for_results_table/all.csv')

#compile()


def create_table(results_file, out_file):
    # Import(s)
    import numpy as np
    import pandas as pd

    # Action
    df = pd.read_csv(results_file)

    raw_ids = np.array(df["ID"])
    ids = np.array([])
    
    for id in raw_ids: 
        id = id.replace('_',' ')
        ids = np.append(ids, id)
    
    ras = np.array(df["RA"])
    decs = np.array(df["DEC"])
    rs = np.array(df["mean_r"])
    classes = np.array([str(i).upper() for i in df["class"]])
    nus = np.array(df["NU"])
    raw_periods = np.array(df["PER1"])
    periods = np.array([])
    for period, qm_class in zip(raw_periods, classes):
        if qm_class == 'P':
            periods = np.append(periods, round(period,3))
        else: 
            periods = np.append(periods, '-')
    qs = np.array(df["Q"])
    ms = np.array(df["M"])
    
    with open(out_file,'a') as f:
        for id, ra, dec, r, qm_class, nu, per, q, m in zip(ids, ras, 
        decs, rs, classes, nus, periods, qs, ms):
            f.write(str(id)+' & '+str(round(ra,5))+' & '+str(round(dec,5))
            +' & '+str(round(float(r),3))+' & '+qm_class+' & '+
            str(round(nu,3))+' & '+per+' & '+str(round(q,3))+' & '+
            str(round(m,3))+r' \\'+'\n')
            
            
    
    
    
results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/Modifying_LaTeX_Code/stuff_for_results_table/all.csv'
out_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/mu/Modifying_LaTeX_Code/stuff_for_results_table/formatted.txt'
create_table(results_file, out_file)
    