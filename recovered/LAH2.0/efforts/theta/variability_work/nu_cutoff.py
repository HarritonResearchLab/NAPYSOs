def create_mags_vs_var_plot(percentile,key,data_path,ref_folder_path,copy_path,report_file):
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    from personalastropy.ysospy.variability import sokolovskyNu
    from scipy import interpolate
    from matplotlib.ticker import AutoMinorLocator
    import os, shutil

    # Action
    ids = list(pd.read_csv(key)['ID'])
    mean_mags = []
    nus = []
    
    for id in ids: 
        data_file = pd.read_csv(data_path.replace('+++',id))
        mags = np.array(data_file['mag'])
        magerrs = np.array(data_file['magerr'])
        mean_mag = np.mean(mags)
        nu = sokolovskyNu(x=mags,xerr=magerrs)
        mean_mags.append(mean_mag)
        nus.append(nu)
    
    mean_mags = np.array(mean_mags)
    nus = np.array(nus)

    min_mean = int(np.floor(min(mean_mags)))
    max_mean = int(np.ceil(max(mean_mags)))

    lower_bounds = list(range(min_mean,max_mean))

    # Quick inizialization of plot
    plt.rcParams['font.family'] = 'serif'
    plt.scatter(mean_mags,nus,s=3,color='cornflowerblue',label='NAP YSOs',marker='o')    
    plt.xlabel('R (mag)')
    plt.ylabel('Peak-to-peak variability ('+r'$\nu$'+')')

    median_bounds = []
    cutoffs = []

    for lower_bound in lower_bounds: 
        upper_bound = lower_bound + 1.0
        mask = np.logical_and(mean_mags<upper_bound,mean_mags>=lower_bound)
        masked_nus = nus[mask]
        cutoff = float(np.percentile(masked_nus,q=[percentile]))
        cutoffs.append(cutoff)
        
        med_bound = np.median([lower_bound,upper_bound])
        median_bounds.append(med_bound)

        # Plot 
        # plt.plot([lower_bound,upper_bound],2*[cutoff],color='blue',lw=1,ls='--')    
    
    
    
    f = interpolate.interp1d(median_bounds,cutoffs,kind='linear',fill_value='extrapolate')
    inter_x = np.linspace(min(mean_mags),max(mean_mags),50)
    inter_y = f(inter_x)
    
    
    list_of_names = os.listdir(ref_folder_path)
    
    ids_below_cutoff = []
    nus_below_cutoff = []
    mean_mags_below_cutoff = []

    for id, mean_mag, nu in zip(ids,list(mean_mags),list(nus)):
        if nu < f(mean_mag):
            for name in list_of_names:
                name_list = name.split(':')
                if id == name_list[2].replace('.png',''):
                    new_loc = copy_path.replace('+++',name)
                    shutil.copyfile(os.path.join(ref_folder_path,name),new_loc)
                    ids_below_cutoff.append(id)
                    mean_mags_below_cutoff.append(mean_mag)
                    nus_below_cutoff.append(nu)
                    break

    zipped_list = list(zip(ids_below_cutoff,nus_below_cutoff,mean_mags_below_cutoff))
    df = pd.DataFrame(data=zipped_list,columns=['ID','NU','MEAN_MAG'])
    df.to_csv(report_file,index=False)

    def finish_plot():
        plt.plot(inter_x,inter_y,color='black',lw=1,ls='--',label=str(percentile)+'th Percentile')
        plt.yscale('log')
        plt.xlim(min(mean_mags),max(mean_mags))
        plt.gca().invert_xaxis()
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
        # Some legend stuff
        leg = plt.legend()
        plt.legend(loc='lower left',fontsize=7,fancybox=False,edgecolor='black',shadow=False)
        for lh in leg.legendHandles: 
            lh.set_alpha(0.4)
        leg.get_frame().set_edgecolor('black')
        plt.gca().tick_params(axis='both',which='both',direction='in')
        plt.show()

    finish_plot()
  

dp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++_r.csv'
key_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/key.csv'

ref_folder = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/eta/q_work_with_owen/actual_q_work/q_diagnostics/plots'
copy = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/variability_work/15_percent/below/+++'
report_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/variability_work/15_percent/below_var_cutoff.csv'
create_mags_vs_var_plot(percentile=15,key=key_path,data_path=dp,ref_folder_path=ref_folder,copy_path=copy,report_file=report_file)


def remove_dim_objects(key,data_path,max_mean_mag,report_file):
    '''
    Remove objects from sample if 
    mean g mag > max_mean_mag.
    '''
    
    # Import(s)
    import pandas as pd
    import numpy as np

    # Action
    key_df = pd.read_csv(key)
    ids = list(key_df['ID'])
    bad_ids = []
    mean_mags = []
    for id in ids: 
        data_file = data_path.replace('+++',id)
        data_df = pd.read_csv(data_file)
        mean_mag = np.mean(np.array(data_df['mag']))
        if mean_mag > max_mean_mag:
            print(id)
            bad_ids.append(id)
            mean_mags.append(mean_mag)

    zipped_list = list(zip(bad_ids,mean_mags))
    df = pd.DataFrame(data=zipped_list,columns=['ID','MEAN_MAG'])
    df.to_csv(report_file,index=False)

def create_final_list(key_file,dim_objects_file,non_variables_file,report_file):
    '''
    Create final list of variable objects that meet var. and 
    brightness criteria
    '''

    # Import(s)
    import pandas as pd

    # Action
    key_df = pd.read_csv(key_file)
    key_ids = list(key_df['ID'])
    dim_df = pd.read_csv(dim_objects_file)
    dim_ids = list(dim_df['ID'])
    non_var_df = pd.read_csv(non_variables_file)
    non_var_ids = list(non_var_df['ID'])


    good_ids = []

    for id in key_ids:
        if id not in dim_ids:
            if id not in non_var_ids:
                good_ids.append(id)

    report_df = pd.DataFrame(data=list(zip(good_ids)),columns=['ID'])
    report_df.to_csv(report_file,index=False)

'''
remove_dim_objects(key=key_path,data_path=dp.replace('r.csv','g.csv'),max_mean_mag=20.8,
                report_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/remove_dim_objects/20.8_mag_cutoff.csv')

create_final_list(key_file=key_path,
                dim_objects_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/remove_dim_objects/20.8_mag_cutoff.csv',
                non_variables_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/variability_work/15_percent/below_var_cutoff.csv',
                report_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/theta/good_ids.csv')
'''