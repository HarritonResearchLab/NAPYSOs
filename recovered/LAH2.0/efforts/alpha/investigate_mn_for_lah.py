def lomb_scarg():
    #Import(s)
    from personalastropy.ysospy.variability import lombScargle
    from astropy.io import ascii
    import numpy as np
    import matplotlib.pyplot as plt

    #Action

    red_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/FHK_286_r.tbl'
    r = ascii.read(red_file,format='ipac',delimiter='|') 
    o_red_dates = list(r['hjd'])
    o_red_mags = list(r['mag'])


    red_dates = []
    for item in o_red_dates:
        item = float(item)
        red_dates.append(item)

    red_mags = []
    for item in o_red_mags:
        item = float(item)
        red_mags.append(item)


    n_r_dates = []
    n_r_mags = []
    for date in red_dates:
        if date < 2458448 or date > 2458456:
            n_r_dates.append(date)
            date_index = red_dates.index(date)
            n_r_mags.append(red_mags[date_index])

    #Plot
    plt.scatter(n_r_dates,n_r_mags,color='indianred',marker='.')
    plt.xlabel('Time (MJD)')
    plt.ylabel('Mag')
    plt.gca().invert_yaxis()
    plt.show()


    lombScargle(id='FHK_176',x=n_r_dates,y=n_r_mags,min_period=2,max_period=400,out_type='FHK_286_ls.png')

def investigate_cluster_effect(r_file,g_file):
    #Import(s)
    from astropy.io import ascii
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import odr
    from personalastropy.ysospy.interpolation import cmdPrep2, returnGoodRegions
    from personalastropy.ysospy.handy_scripts import sortData

    #Action
    r = ascii.read(r_file,format='ipac',delimiter='|')
    r_dates = list(r['mjd'])
    r_mags = list(r['mag'])
    r_magerrs = list(r['magerr'])

    g = ascii.read(g_file,format='ipac',delimiter='|')
    g_dates = list(g['mjd'])
    g_mags = list(g['mag'])
    g_magerrs = list(g['magerr'])

    #Get rid of clusters 

    c_r_dates = [] #c for cleaned
    c_r_mags = []
    c_r_magerrs = []

    for item in r_dates: 
        if item > 58452 or item < 58448:
            item_index = r_dates.index(item)
            c_r_dates.append(item)
            c_r_mags.append(r_mags[item_index])
            c_r_magerrs.append(r_magerrs[item_index])

    #CMD stuff
    red_sorter = sortData(x=c_r_dates,y=[c_r_mags,c_r_magerrs])
    srd = red_sorter[0]
    srm = red_sorter[1][0]
    sre = red_sorter[1][1]

    green_sorter = sortData(x=g_dates,y=[g_mags,g_magerrs])

    sgd = green_sorter[0]
    sgm = green_sorter[1][0]
    sge = green_sorter[1][1] 

    gis = returnGoodRegions(x=sgd,y=sgm,max_sep=3,min_card=2) #good green intervals func
    
    green_date_intervals = gis[1]

    cmdData = cmdPrep2(x=sgd,y=sgm,magerrs=sge,gdis=green_date_intervals,x2=srd,y2=srm,magerrs2=sre,kind='linear')

    cmd_dates = list(cmdData[0])
    cmd_x = list(cmdData[2])
    cmd_y = list(cmdData[1])
    cmd_xerr = list(cmdData[3])
    cmd_yerr = list(cmdData[4])

    #Plot
    plt.rcParams['font.family'] = 'serif'
    fig, axs = plt.subplots(2,3)
    

    axs[0,0].errorbar(srd,srm,yerr=sre,color='tomato',marker='.',lw=0,elinewidth=0.3,ms=1,label='640 nm')
    axs[0,0].errorbar(sgd,sgm,yerr=sge,color='mediumseagreen',marker='.',lw=0,ms=1,elinewidth=0.3,label='480 nm')
    axs[0,0].set_xlabel('time (MJD)')
    axs[0,0].set_ylabel('mag')
    axs[0,0].invert_yaxis()
    box = axs[0,0].get_position()
    axs[0,0].set_position([box.x0,box.y0,box.width,box.height*0.7])
    axs[0,0].legend(bbox_to_anchor=(0.5,1.1),loc='center',fontsize='xx-small')

    axs[0,1].errorbar(cmd_dates,cmd_x,yerr=cmd_xerr,color='cornflowerblue',marker='.',lw=0,elinewidth=0.3,ms=1)
    axs[0,1].set_xlabel('time (MJD)')
    axs[0,1].set_ylabel('[480]-[640]')
    axs[0,1].invert_yaxis()
    

    axs[0,2].errorbar(cmd_x,cmd_y,xerr=cmd_xerr,yerr=cmd_yerr,ms=0,elinewidth=0.5,alpha=0.4,lw=0,color='lightgrey')
    axs[0,2].scatter(cmd_x,cmd_y,color='cornflowerblue',marker='.')
    
    #ir slope setup 
    irs = 3.81
    max_y = max(cmd_y)
    min_y = min(cmd_y)
    med_y = float(np.median(cmd_y))
    med_x = float(np.median(cmd_x))

    max_x = ((max_y-med_y)/(irs))+med_x
    min_x = ((min_y-med_y)/(irs))+med_x

    x_vals = np.array([min_x,max_x,med_x])
    y_vals = (irs*(x_vals-med_x))+med_y

    axs[0,2].plot(x_vals,y_vals,color='tomato',label='IR slope: 3.81',ls=':')

    
    #exclude < 2.25th and > 97.5th percentile values


    g_percentiles = np.percentile(cmd_y,q=[2.5,97.5])
    g_minus_r_percentiles = np.percentile(cmd_x,[2.5,97.5])

    for item in cmd_y:
        if item < g_percentiles[0] or item > g_percentiles[1]:
            item_index = cmd_y.index(item)
            cmd_y.remove(item)
            cmd_x.remove(cmd_x[item_index])
            cmd_xerr.remove(cmd_xerr[item_index])
            cmd_yerr.remove(cmd_yerr[item_index])

    for item in cmd_x:
        if item < g_minus_r_percentiles[0] or item > g_minus_r_percentiles[1]:
            item_index = cmd_x.index(item)
            cmd_x.remove(item)
            cmd_y.remove(cmd_y[item_index])
            cmd_xerr.remove(cmd_xerr[item_index])
            cmd_yerr.remove(cmd_yerr[item_index])
    
    g_vals = np.array(cmd_y)
    g_minus_r_magerrs=np.array(cmd_xerr)
    cmd_g_magerrs=np.array(cmd_yerr)
    g_minus_r = np.array(cmd_x)

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

    axs[1,0].errorbar(r_dates,r_mags,yerr=r_magerrs,color='tomato',marker='.',lw=0,elinewidth=0.3,ms=1,label='640 nm')
    #axs[1,0].errorbar(sgd,sgm,yerr=sge,color='mediumseagreen',marker='.',lw=0,ms=1,elinewidth=0.3,label='480 nm')
    axs[1,0].set_xlabel('time (MJD)')
    axs[1,0].set_ylabel('mag')
    axs[1,0].invert_yaxis()
    box = axs[1,0].get_position()
    axs[1,0].set_position([box.x0,box.y0,box.width,box.height*0.7])
    axs[1,0].legend(bbox_to_anchor=(0.5,1.1),loc='center',fontsize='xx-small')
    


    axs[1,1].set_axis_off()
    axs[1,2].set_axis_off()
    plt.subplots_adjust(wspace=0.6,hspace=0.5)
    plt.show()


investigate_cluster_effect(r_file='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/2MASS_J20505568+4421442_r.tbl',g_file='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/2MASS_J20505568+4421442_g.tbl')


