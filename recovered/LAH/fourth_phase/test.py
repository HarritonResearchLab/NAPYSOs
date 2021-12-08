#imports
from personalastropy.ysospy.interpolation import returnGoodRegions, cmdPrep
from scipy import interpolate
from personalastropy.ysospy.handy_scripts import sortData
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import seaborn as sns

def cmd_stuff():
    #Import(s)
    from astropy.io import ascii
    #Action

    coordinates = []
    with open('/home/thaddaeus/FMU/HRL/LAH/third_phase/data/obsids.txt','r') as f:
        for line in f:
            line = line.replace('\n','')
            coordinates.append(line)

    for coordinate in coordinates: 
        path_temp = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_+++++_=.tbl'
        red_file = path_temp.replace('+++++',coordinate)
        green_file = red_file.replace('=','g')
        red_file=red_file.replace('=','r')

        r=ascii.read(red_file,format='ipac',delimiter='|')
        g=ascii.read(green_file,format='ipac',delimiter='|')

        red_mags=r['mag']
        red_dates=r['hjd']
        red_magerrs=r['magerr']

        green_mags = g['mag']
        green_dates=g['hjd']
        green_magerrs=g['magerr']


        if len(red_mags) > 5 and len(green_mags) > 5:
            func = sortData(x=red_dates,y=[red_mags,red_magerrs])
            srd = func[0]
            srm=func[1][0]
            sre=func[1][1]
            gfunc=sortData(x=green_dates,y=[green_mags,green_magerrs])
            sgd = gfunc[0]
            sgm=gfunc[1][0]
            sge=gfunc[1][1]
            gis = returnGoodRegions(x=sgd,y=sgm,max_sep=3,min_card=2)
            green_date_intervals = gis[1]
            green_mag_intervals = gis[2]


            interf = cmdPrep(sgd,sgm,green_date_intervals,srd,srm,'linear')
            
            #do the lin regress
            r_func = linregress(interf[2],interf[1])
            
            slope = r_func[0]
            intercept = r_func[1]
            r=r_func[2]
            r_sq=r**2
            p=r_func[3]

            if r_sq > 0.75:
                print(coordinate)
                print(r_sq)

            reg_linex=np.linspace(min(interf[2]),max(interf[2]),5)
            
            plt.scatter(interf[2],interf[1],color='cornflowerblue',s=3,marker='o')
            sns.regplot(x=interf[2],y=interf[1],label=(r'$r^2$'+': '+str(round(r_sq,2))+'\nslope: ')+str(round(slope,2)),fit_reg=True,scatter=False,ci=95,order=1)
            plt.gca().invert_yaxis()
            plt.xlabel('g-r')
            plt.ylabel('g')
            plt.title('Color-Mag Diagram')
            plt.legend()

            
            
            #plt.subplots_adjust(hspace=0.4)
            #plt.show()
            
            plt.savefig('/home/thaddaeus/FMU/HRL/LAH/fourth_phase/color_mag_work/lcmg_with_reg/'+coordinate+'.png')
            
            plt.clf()

def comparisons():
    #Import(s)
    from personalastropy.ysospy.variability import sokolovskyNu
    import scipy.stats as stats
    import seaborn as sns
    from astropy.io import ascii

    #Action

    coordinates = []
    with open('/home/thaddaeus/FMU/HRL/LAH/third_phase/data/obsids.txt','r') as f:
        for line in f:
            coordinates.append(line.replace('\n',''))
    
    
    skews = []
    kurts = []
    nus = []
    stdevs = []

    for coordinate in coordinates:
        red_file = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_++++++++++_r.tbl'.replace('++++++++++',coordinate)
        r=ascii.read(red_file,format='ipac',delimiter='|')

        red_mags=np.array(r['mag'])
        red_magerrs=np.array(r['magerr'])

        if len(red_mags) > 3 and len(red_magerrs) > 3: 
            skews.append(stats.skew(red_mags))
            kurts.append(stats.kurtosis(red_mags))
            stdevs.append(np.std(red_mags,ddof=1))
            nus.append(sokolovskyNu(red_mags,red_magerrs))

    sns.set_style('darkgrid')
    
    fig, axs = plt.subplots(2, 3)
    

    axs[0,0].scatter(skews,kurts,s=3)
    axs[0,0].set_xlabel('Skew')
    axs[0,0].set_ylabel('Kurtosis')

    axs[0,1].scatter(skews,stdevs,s=3)
    axs[0,1].set_xlabel('Skew')
    axs[0,1].set_ylabel(r'$\sigma$')

    axs[0,2].scatter(skews,nus,s=3)
    axs[0,2].set_xlabel('Skew')
    axs[0,2].set_ylabel('Sokolovsky '+r'$\nu$')

    axs[1,0].scatter(stdevs,kurts,s=3)
    axs[1,0].set_xlabel(r'$\sigma$')
    axs[1,0].set_ylabel('Kurtosis')

    axs[1,1].scatter(nus,kurts,s=3)
    axs[1,1].set_xlabel('Sokolovsky '+r'$\nu$')
    axs[1,1].set_ylabel('Kurtosis')

    axs[1,2].scatter(nus,stdevs,s=3)
    axs[1,2].set_xlabel('Sokolovsky '+r'$\nu$')
    axs[1,2].set_ylabel(r'$\sigma$')

    plt.subplots_adjust(hspace=0.4,wspace=0.8)
    plt.show()


cmd_stuff()

    
