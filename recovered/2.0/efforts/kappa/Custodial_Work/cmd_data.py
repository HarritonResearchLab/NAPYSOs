def create_cmd_data(key,r_path,g_path,cmd_csv_temp): 
    #Import(s)
    import pandas as pd
    from personalastropy.ysospy.interpolation import cmdPrep2, returnGoodRegions
    from personalastropy.ysospy.handy_scripts import sortData
    from progress.bar import Bar 

    #Action

    df = pd.read_csv(key)

    ids = list(df['ID'])
    bar = Bar('Processing',max=len(ids))
    for id in ids:
        r_file = r_path.replace('+++',id)
        r = pd.read_csv(r_file)
        r_mjds = list(r['mjd'])
        r_mags = list(r['mag'])
        r_magerrs = list(r['magerr'])

        g_file = g_path.replace('+++',id)
        g = pd.read_csv(g_file)
        g_mjds = list(g['mjd'])
        g_mags = list(g['mag'])
        g_magerrs = list(g['magerr'])

        red_sorter = sortData(x=r_mjds,y=[r_mags,r_magerrs])

        srd = red_sorter[0]
        srm = red_sorter[1][0]
        sre = red_sorter[1][1]

        green_sorter = sortData(x=g_mjds,y=[g_mags,g_magerrs])

        sgd = green_sorter[0]
        sgm = green_sorter[1][0]
        sge = green_sorter[1][1] 

        #Process and save cmd data to csv 
        gis = returnGoodRegions(x=sgd,y=sgm,max_sep=3,min_card=2) #good green intervals func
        
        green_date_intervals = gis[1]
        
        cmdData = cmdPrep2(x=sgd,y=sgm,magerrs=sge,gdis=green_date_intervals,x2=srd,y2=srm,magerrs2=sre,kind='linear')
        
        cmd_csv_file = cmd_csv_temp.replace('+++',id)
        
        date_col = list(cmdData[0])
        col_1 = list(cmdData[2])
        col_2 = list(cmdData[1])
        col_3 = list(cmdData[3])
        col_4 = list(cmdData[4])
        col_5 = list(cmdData[5]) #red dates
        col_6 = list(cmdData[6]) #red mags
        col_7 = list(cmdData[7]) #red magerrs
    
        print(cmdData[5])
        
        cmd_zipped_list = list(zip(date_col,col_1,col_2,col_3,col_4,col_5,col_6,col_7))

        cmd_df = pd.DataFrame(cmd_zipped_list,columns=['g_date','g-r_val','g_mag','g-r_magerr','g_magerr','r_date','r_mag','r_magerr'])   

        cmd_df.to_csv(cmd_csv_file,index=False)
        
        bar.next()

    bar.finish()

### RUN IT

key_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/key.csv'
r_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
g_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_g.csv'
cmd_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/cmd/+++.csv'
create_cmd_data(key=key_file,r_path=r_path,g_path=g_path,cmd_csv_temp=cmd_temp)

'''

# ONLY FOR THE FHK_176 RECOVERY!!!!! 
# !!!! 
# !!!!
key_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/fake_key.csv'
r_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/light_curves/+++_r.csv'
g_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/light_curves/+++_g.csv'
cmd_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/cmds/+++.csv'
#create_cmd_data(key=key_file,r_path=r_path,g_path=g_path,cmd_csv_temp=cmd_temp)


def make_missing_cmd(old_cmd_data, new_r_lc, save_dir):
    # Import(s)
    import numpy as np
    import pandas as pd
    
    # Action
    cmd_dates = np.array(pd.read_csv(old_cmd_data)['date'])
    lc_dates = np.array(pd.read_csv(new_r_lc)['mjd'])
    merged_dates = np.intersect1d(cmd_dates, lc_dates)
    
    fake_df = pd.DataFrame(list(zip(lc_dates)), columns=['date'])
    
    cmd_df = pd.read_csv(old_cmd_data)
    new_cmd = cmd_df.merge(fake_df, on='date')
    new_cmd.to_csv(save_dir+'/FHK_176.csv', index=False)
        
old_cmd_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/cmd_data/FHK_176.csv'
make_missing_cmd(old_cmd_file, r_path.replace('+++', 'FHK_176'), 
'/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/')
'''