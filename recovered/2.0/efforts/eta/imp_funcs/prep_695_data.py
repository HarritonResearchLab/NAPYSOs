def clean_clusters(dates,paired_lists,tolerance):
    #Import(s)
    import numpy as np
    import math

    #Action

    min_date = math.floor(np.min(dates))
    max_date = math.ceil(np.max(dates))
    date_bins = list(range(min_date,max_date,1))
    binned_dates = np.digitize(x=dates,bins=date_bins)
    bin_counts = np.bincount(binned_dates)

    nz_bin_counts = bin_counts[np.nonzero(bin_counts)] #Get rid of bin_counts = 0 for calculations

    mean_counts = np.mean(nz_bin_counts)
    std_counts = np.std(nz_bin_counts,ddof=1)

    cutoff = mean_counts+float(tolerance)*(std_counts)

    intervals_to_remove = []
    
    bin_counts = list(bin_counts)

    for count in bin_counts:
        if count > cutoff:
            bin_index = bin_counts.index(count)
            lower_edge = date_bins[bin_index-1]
            upper_edge = date_bins[bin_index]
            edge_string = str(lower_edge)+':'+str(upper_edge)
            intervals_to_remove.append(edge_string)

    indices_to_remove = []
    intervals_to_remove = list(set(intervals_to_remove))
    for interval in intervals_to_remove:
        
        interval_list = interval.split(':')
        lower = float(interval_list[0])
        upper = float(interval_list[1])
        for date in dates:
            if date > lower and date < upper:
                date_index = dates.index(date)
                indices_to_remove.append(date_index)

    for i in sorted(indices_to_remove, reverse=True):
        del dates[i]
        for paired_list in paired_lists:
            del paired_list[i]

    return dates, paired_lists

def rename_red_files(folder):
    #Import(s)
    import os

    #Action
    
    for file in os.listdir(folder):
        if file[-4:] == '.tbl':
            old_file = os.path.join(folder,file)    
            new_file = old_file.replace('.tbl','_r.tbl')
            os.rename(old_file,new_file)
    
def process_and_plot(key,r_path,g_path,placeholder,cmd_csv_dir,g_csv_dir,r_csv_dir,plot_dir): ###THIS FUNC HAS BEEN UPDATED TO RETURN RED ARRAY COLUMN AS WELL
    #Import(s)
    import pandas as pd
    import numpy as np
    from astropy.io import ascii
    from personalastropy.ysospy.interpolation import cmdPrep2, returnGoodRegions
    from personalastropy.ysospy.handy_scripts import sortData
    import matplotlib.pyplot as plt
    from progress.bar import Bar 

    #Action

    df = pd.read_csv(key)

    ids = list(df['ID'])
    bar = Bar('Processing',max=len(ids))
    for id in ids:
        r_file = r_path.replace(placeholder,id)
        r = ascii.read(r_file,format='ipac',delimiter='|')
        r_mjds = list(r['mjd'])
        r_mags = list(r['mag'])
        r_magerrs = list(r['magerr'])

        g_file = g_path.replace(placeholder,id)
        g = ascii.read(g_file,format='ipac',delimiter='|')
        g_mjds = list(g['mjd'])
        g_mags = list(g['mag'])
        g_magerrs = list(g['magerr'])


        cleaned_r_mjds = []
        cleaned_r_mags = []
        cleaned_r_magerrs = []

        for item in r_mjds: 
            if item > 58456 or item < 58448:
                item_index = r_mjds.index(item)
                cleaned_r_mjds.append(item)
                cleaned_r_mags.append(r_mags[item_index])
                cleaned_r_magerrs.append(r_magerrs[item_index])

        red_sorter = sortData(x=cleaned_r_mjds,y=[cleaned_r_mags,cleaned_r_magerrs])

        srd = red_sorter[0]
        srm = red_sorter[1][0]
        sre = red_sorter[1][1]
        
        '''
        cleaned_green = clean_clusters(dates=g_mjds,paired_lists=[g_mags,g_magerrs],tolerance=4)
        cleaned_g_mjds = cleaned_green[0]
        cleaned_g_mags = cleaned_green[1][0]
        cleaned_g_magerrs = cleaned_green[1][1]
        green_sorter = sortData(x=cleaned_g_mjds,y=[cleaned_g_mags,cleaned_g_magerrs])
        '''

        green_sorter = sortData(x=g_mjds,y=[g_mags,g_magerrs])

        sgd = green_sorter[0]
        sgm = green_sorter[1][0]
        sge = green_sorter[1][1] 

        '''
        #Save cleaned red and green data to csv files
        
        g_zipped_list = list(zip(sgd,sgm,sge))
        g_csv_file = g_csv_dir.replace(placeholder,id)
        g_df = pd.DataFrame(g_zipped_list,columns=['mjd','mag','magerr'])
        g_df.to_csv(g_csv_file,index=False)

        r_zipped_list = list(zip(srd,srm,sre))
        r_csv_file = r_csv_dir.replace(placeholder,id)
        r_df = pd.DataFrame(r_zipped_list,columns=['mjd','mag','magerr'])
        r_df.to_csv(r_csv_file,index=False)
        '''
        #Process and save cmd data to csv 
        gis = returnGoodRegions(x=sgd,y=sgm,max_sep=3,min_card=2) #good green intervals func
        
        green_date_intervals = gis[1]

        cmdData = cmdPrep2(x=sgd,y=sgm,magerrs=sge,gdis=green_date_intervals,x2=srd,y2=srm,magerrs2=sre,kind='linear')

        cmd_csv_file = cmd_csv_dir.replace(placeholder,id)
        
        date_col = list(cmdData[0])
        col_1 = list(cmdData[2])
        col_2 = list(cmdData[1])
        col_3 = list(cmdData[3])
        col_4 = list(cmdData[4])
        col_5 = list(cmdData[5]) #red dates
        col_6 = list(cmdData[6]) #red mags
        col_7 = list(cmdData[7]) #red magerrs

        
        cmd_zipped_list = list(zip(date_col,col_1,col_2,col_3,col_4,col_5,col_6,col_7))

        cmd_df = pd.DataFrame(cmd_zipped_list,columns=['g_date','g-r_val','g_mag','g-r_magerr','g_magerr','r_date','r_mag','r_magerr'])   

        cmd_df.to_csv(cmd_csv_file,index=False)
        
        def Plot_lightcurve_and_cmd(): 
            '''
            fig, axs = plt.subplots(2,1)

            #Lightcurve
            sns.set_style('white')
            axs[0].errorbar(srd,srm,yerr=sre,color='tomato',marker='.',lw=0,elinewidth=0.3,ms=1,label='640 nm')
            axs[0].errorbar(sgd,sgm,yerr=sge,color='mediumseagreen',marker='.',lw=0,ms=1,elinewidth=0.3,label='480 nm')
            axs[0].set_xlabel('time (MJD)')
            axs[0].set_ylabel('mag')
            axs[0].invert_yaxis()
            box = axs[0].get_position()
            axs[0].set_position([box.x0,box.y0,box.width*0.7,box.height])
            axs[0].legend(bbox_to_anchor=(1,0.5),loc='center left',fontsize='xx-small')

            #CMD

            #axs[1].errorbar(col_1,col_2,xerr=col_3,yerr=col_4,lw=0,elinewidth=0,ecolor='darkgrey',marker='.',mfc='cornflowerblue')
            irs = 3.81
            axs[1].scatter(col_1,col_2,color='cornflowerblue',marker='.')
            
            #Plot the interstellar reddening line
            max_y = max(col_2)
            min_y = min(col_2)
            med_y = float(np.median(col_2))
            med_x = float(np.median(col_1))

            max_x = ((max_y-med_y)/(irs))+med_x
            min_x = ((min_y-med_y)/(irs))+med_x

            x_vals = np.array([min_x,max_x,med_x])
            y_vals = (irs*(x_vals-med_x))+med_y

            axs[1].plot(x_vals,y_vals,color='indianred',label='IR slope',ls='--')
            
            #Finish up axs[1]
            axs[1].set_xlabel('[480] - [640]')
            axs[1].set_ylabel('[480]')
            axs[1].invert_yaxis()
            box = axs[1].get_position()
            axs[1].set_position([box.x0,box.y0,box.width*0.7,box.height])
            axs[1].legend(bbox_to_anchor=(1,0.5),loc='center left',fontsize='xx-small')

            plt.subplots_adjust(hspace=0.3)

            plt.savefig(plot_dir.replace(placeholder,id))
            fig.clf()
            plt.close()
            '''
        
        bar.next()

def get_rid_off_bad_files(key,csv_file):
    #Import(s)
    import pandas as pd
    import shutil
    from astropy.io import ascii

    #Action

    df = pd.read_csv(key)

    ids = list(df['ID'])
    ras = list(df['RA'])
    decs = list(df['DEC'])

    new_ids = []
    new_ras = []
    new_decs = []

    for id in ids:
        r_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/ZTF_R/+++_r.tbl'.replace('+++',id)
        r = ascii.read(r_file,format='ipac',delimiter='|')
        r_mjds = list(r['mjd'])

        g_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/ZTF_G/+++_g.tbl'.replace('+++',id)
        g = ascii.read(g_file,format='ipac',delimiter='|')
        g_mjds = list(g['mjd'])

        if len(r_mjds) > 20 and len(g_mjds) > 20:
            id_index = ids.index(id)
            new_ids.append(id)
            new_ras.append(ras[id_index])
            new_decs.append(decs[id_index])

    zipped_list = list(zip(new_ids,new_ras,new_decs))

    df = pd.DataFrame(zipped_list,columns=['ID','RA','DEC'])    
    df.to_csv(csv_file,index=False)

    df = pd.read_csv(csv_file)

    ids = list(df['ID'])
    
    for id in ids: 
        old_r = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/ZTF_R/+++_r.tbl'.replace('+++',id)
        old_g = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/ZTF_G/+++_g.tbl'.replace('+++',id)  
        new_r = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++_r.tbl'.replace('+++',id)
        new_g = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++_g.tbl'.replace('+++',id)
        shutil.copyfile(old_r,new_r)
        shutil.copyfile(old_g,new_g)

def ols_on_cmds(key,cmd_path,g_path,r_path,placeholder,out_file):
    #Import(s)
    from scipy import stats
    import pandas as pd
    import numpy as np
    from astropy.io import ascii
    from personalastropy.ysospy.variability import sokolovskyNu
    from progress.bar import Bar

    #Action
    df = pd.read_csv(key)
    ids = list(df['ID'])

    good_ids = []
    r_stds = []
    g_stds = []
    r_nus = []
    g_nus = []
    slopes = []
    r_sqs = []
    
    bar = Bar('Processing',max=len(ids))
    for id in ids:
        #Lightcurve statistics
        r_file = r_path.replace(placeholder,id)
        g_file = g_path.replace(placeholder,id)

        g = ascii.read(g_file,format='ipac',delimiter='|')
        r = ascii.read(r_file,format='ipac',delimiter='|')

        g_mags = np.array(g['mag'])
        g_magerrs = np.array(g['magerr'])
        r_mags = np.array(r['mag'])
        r_magerrs = np.array(r['magerr'])

        if len(r_mags) >=30 and len(g_mags) >=30:
            
            good_ids.append(id)
            r_stds.append(np.std(r_mags,ddof=1))
            g_stds.append(np.std(g_mags,ddof=1))

            r_nus.append(sokolovskyNu(r_mags,r_magerrs))
            g_nus.append(sokolovskyNu(g_mags,g_magerrs))

            #CMD regression
            cmd_file = cmd_path.replace(placeholder,id)
            cmd_df = pd.read_csv(cmd_file)

            x = cmd_df['g-r']
            y = cmd_df['g']

            ols = stats.linregress(x,y)
            slope = ols[0]
            r_sq = (ols[2])**2

            slopes.append(slope)
            r_sqs.append(r_sq)

        bar.next()
    
    zipped_list = list(zip(good_ids,slopes,r_sqs,r_stds,r_nus,g_stds,g_nus))

    df_to_save = pd.DataFrame(zipped_list,columns=['ID','Slope','R^2','r std','r nu','g std','g nu'])   

    df_to_save.to_csv(out_file,index=False)

def analyze_ols(ols_results):
    #Import(s)
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    #Actions
    df = pd.read_csv(ols_results)
    slopes = list(df['Slope'])
    r_sqs = list(df['R^2'])
    r_stds = list(df['r std'])
    r_nus = list(df['r nu'])
    g_stds = list(df['g std'])
    g_nus = list(df['g nu'])

    m80_slopes = [] #r^2 > 0.8
    m80_rsqs = []
    r80_stds = []
    r80_nus = []
    g80_stds = []
    g80_nus = []

    #get r_sq > 80 values
    for item in r_sqs: 
        if item > 0.8:
            item_index = r_sqs.index(item)
            m80_rsqs.append(r_sqs[item_index])
            m80_slopes.append(slopes[item_index])
            r80_stds.append(r_stds[item_index])
            r80_nus.append(r_nus[item_index])
            g80_stds.append(g_stds[item_index])
            g80_nus.append(g_nus[item_index])

    #Create plot
    fig, axs = plt.subplots(2,3)
    sns.set_style('white')
    
    slope_label = 'Slope ('+r'$r^2>0.8$'+')'

    axs[0,0].hist(m80_slopes,bins=18,color='cornflowerblue')
    axs[0,0].set_ylabel('Frequency')    
    axs[0,0].set_xlabel(slope_label)
 
    axs[1,0].scatter(r_sqs,slopes,color='cornflowerblue',marker='.')
    axs[1,0].set_ylabel('Slope')
    axs[1,0].set_xlabel(r'$r^2$')

    axs[0,1].scatter(m80_rsqs,r80_stds,color='tomato',marker='.',label='[640]')
    axs[0,1].set_ylabel(r'$\sigma$')
    axs[0,1].set_xlabel(slope_label)
    axs[0,1].legend(loc='upper center',bbox_to_anchor=(0.5,1.2),fontsize=6)

    axs[1,1].scatter(m80_rsqs,g80_stds,color='seagreen',marker='.',label='[420]')
    axs[1,1].set_ylabel(r'$\sigma$')
    axs[1,1].set_xlabel(slope_label)
    axs[1,1].legend(loc='upper center',bbox_to_anchor=(0.5,1.2),fontsize=6)

    axs[0,2].scatter(m80_rsqs,r80_nus,color='tomato',marker='.',label='[640]')
    axs[0,2].set_ylabel(r'$\nu$')
    axs[0,2].set_xlabel(slope_label)
    axs[0,2].legend(loc='lower right',fontsize=6)

    axs[1,2].scatter(m80_rsqs,g80_nus,color='seagreen',marker='.',label='[420]')
    axs[1,2].set_ylabel(r'$\nu$')
    axs[1,2].set_xlabel(slope_label)
    axs[1,2].legend(loc='upper center',bbox_to_anchor=(0.5,1.2),fontsize=6)

    plt.subplots_adjust(wspace=0.5,hspace=0.4)
    plt.show()

def demonstrate_for_MN(cmd_data,r_data,g_data):
    #Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from astropy.io import ascii
    from scipy import odr 

    #Action

    df = pd.read_csv(cmd_data)
    cmd_dates = np.array(df['date'])
    g_minus_r = np.array(df['g-r'])
    g_minus_r_magerrs = np.array(df['g-r magerr'])
    g_vals = np.array(df['g'])
    cmd_g_magerrs = np.array(df['g magerr'])

    r = ascii.read(r_data,format='ipac',delimiter='|')
    r_dates = np.array(r['mjd'])
    r_mags = np.array(r['mag'])
    r_magerrs = np.array(r['magerr'])

    g = ascii.read(g_data,format='ipac',delimiter='|')
    g_dates = np.array(g['mjd'])
    g_mags = np.array(g['mag'])
    g_magerrs = np.array(g['magerr'])


    #Plot
    #sns.set_style('white')
    plt.rcParams['font.family'] = 'serif'
    fig, axs = plt.subplots(2,3)
    

    axs[0,0].errorbar(r_dates,r_mags,yerr=r_magerrs,color='tomato',marker='.',lw=0,elinewidth=0.3,ms=1,label='640 nm')
    axs[0,0].errorbar(g_dates,g_mags,yerr=g_magerrs,color='mediumseagreen',marker='.',lw=0,ms=1,elinewidth=0.3,label='480 nm')
    axs[0,0].set_xlabel('time (MJD)')
    axs[0,0].set_ylabel('mag')
    axs[0,0].invert_yaxis()
    box = axs[0,0].get_position()
    axs[0,0].set_position([box.x0,box.y0,box.width,box.height*0.7])
    axs[0,0].legend(bbox_to_anchor=(0.5,1.1),loc='center',fontsize='xx-small')

    axs[0,1].errorbar(cmd_dates,g_minus_r,yerr=g_minus_r_magerrs,color='cornflowerblue',marker='.',lw=0,elinewidth=0.3,ms=1)
    axs[0,1].set_xlabel('time (MJD)')
    axs[0,1].set_ylabel('[480]-[640]')
    axs[0,1].invert_yaxis()
    

    axs[0,2].errorbar(g_minus_r,g_vals,xerr=g_minus_r_magerrs,yerr=cmd_g_magerrs,ms=0,elinewidth=0.5,alpha=0.4,lw=0,color='lightgrey')
    axs[0,2].scatter(g_minus_r,g_vals,color='cornflowerblue',marker='.')
    
    #ir slope setup 
    irs = 3.81
    max_y = max(g_vals)
    min_y = min(g_vals)
    med_y = float(np.median(g_vals))
    med_x = float(np.median(g_minus_r))

    max_x = ((max_y-med_y)/(irs))+med_x
    min_x = ((min_y-med_y)/(irs))+med_x

    x_vals = np.array([min_x,max_x,med_x])
    y_vals = (irs*(x_vals-med_x))+med_y

    axs[0,2].plot(x_vals,y_vals,color='tomato',label='IR slope: 3.81',ls=':')

    
    #exclude < 2.25th and > 97.5th percentile values
    g_vals = list(g_vals)
    g_minus_r_magerrs=list(g_minus_r_magerrs)
    cmd_g_magerrs=list(cmd_g_magerrs)
    g_minus_r = list(g_minus_r)

    g_percentiles = np.percentile(g_vals,q=[2.25,97.5])
    g_minus_r_percentiles = np.percentile(g_minus_r,[2.25,97.5])

    for item in g_vals:
        if item < g_percentiles[0] or item > g_percentiles[1]:
            item_index = g_vals.index(item)
            g_vals.remove(item)
            g_minus_r.remove(g_minus_r[item_index])
            cmd_g_magerrs.remove(cmd_g_magerrs[item_index])
            g_minus_r_magerrs.remove(g_minus_r_magerrs[item_index])

    for item in g_minus_r:
        if item < g_minus_r_percentiles[0] or item > g_minus_r_percentiles[1]:
            item_index = g_minus_r.index(item)
            g_minus_r.remove(item)
            g_vals.remove(g_vals[item_index])
            cmd_g_magerrs.remove(cmd_g_magerrs[item_index])
            g_minus_r_magerrs.remove(g_minus_r_magerrs[item_index])
    
    g_vals = np.array(g_vals)
    g_minus_r_magerrs=np.array(g_minus_r_magerrs)
    cmd_g_magerrs=np.array(g_magerrs)
    g_minus_r = np.array(g_minus_r)

    #odr 
    def line(x,a,b):
        y = (a*x)+b
        return y
    def odr_line(B,x):
        y = B[0]*x+B[1]
        return y
    def perform_odr(x,y,xerr,yerr):
        linear = odr.Model(odr_line)
        mydata = odr.Data(x,y,wd=1./xerr,we=1./yerr)
        myodr = odr.ODR(mydata,linear,beta0=[0,0])
        output = myodr.run()
        return output

    odr_regression = perform_odr(x=g_minus_r,y=g_vals,xerr=g_minus_r_magerrs,yerr=cmd_g_magerrs)
    odr_slope = odr_regression.beta[0]
    odr_slope_error = odr_regression.sd_beta[0]
    error = 100*odr_slope_error/odr_slope #in percent
    ods_error_label = 'Slope: '+str(round(odr_slope,3))+'\nError: '+str(round(error,3))
    ods_simple_label = 'odr Slope: '+str(round(odr_slope,3))
    axs[0,2].plot(g_minus_r,line(g_minus_r,odr_regression.beta[0],odr_regression.beta[1]),color='mediumseagreen',ls=':',label=(ods_simple_label))
    
    #Finish up axs[0,2]
    axs[0,2].set_xlabel('[480] - [640]')
    axs[0,2].set_ylabel('[480]')
    axs[0,2].invert_yaxis()
    box = axs[0,2].get_position()
    axs[0,2].set_position([box.x0,box.y0,box.width,box.height*0.7])
    axs[0,2].legend(bbox_to_anchor=(0.5,1.1),loc='center',fontsize='xx-small')


    axs[1,0].set_axis_off()
    axs[1,1].set_axis_off()
    axs[1,2].set_axis_off()
    plt.subplots_adjust(wspace=0.4)
    plt.show()


    
    
    




r_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++++_r.tbl'
g_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++++_g.tbl'
csv_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/cmd_v2/+++++.csv'
key_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/key.csv'

g_csv_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++++_g.csv'
r_csv_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/data/lightcurves/+++++_r.csv'

process_and_plot(key=key_path,r_path=r_temp,g_path=g_temp,placeholder='+++++',cmd_csv_dir=csv_temp,g_csv_dir='',r_csv_dir='',plot_dir='')