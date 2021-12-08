#Import(s)
import pandas as pd

#Functions

def create_plots(IDs,file_temp,out_temp):
    #Import(s)
    import matplotlib.pyplot as plt
    import pandas as pd
    from astropy.io import ascii
    from progress.bar import Bar

    #Action
    bar = Bar('Plotting',max=len(IDs))

    for id in IDs:

        data_temp = file_temp.replace('+++',id)
        g_file = data_temp.replace('*','g')
        r_file = data_temp.replace('*','r')

        r = ascii.read(r_file,format='ipac',delimiter='|')
        r_dates = list(r['mjd'])
        r_mags = list(r['mag'])
        r_magerrs = list(r['magerr'])

        g = ascii.read(g_file,format='ipac',delimiter='|')
        g_dates = list(g['mjd'])
        g_mags = list(g['mag'])
        g_magerrs = list(g['magerr'])

        c_r_dates = [] #c for cleaned
        c_r_mags = []
        c_r_magerrs = []

        for item in r_dates: 
            if item > 58456 or item < 58448:
                item_index = r_dates.index(item)
                c_r_dates.append(item)
                c_r_mags.append(r_mags[item_index])
                c_r_magerrs.append(r_magerrs[item_index])

        #Plot
        plt.rcParams['font.family'] = 'serif'
        
        fig, axs = plt.subplots(2,1)

        axs[0].errorbar(c_r_dates,c_r_mags,yerr=c_r_magerrs,color='tomato',marker='.',ms=0.5,lw=0,elinewidth=0.5,label='red')
        axs[0].errorbar(g_dates,g_mags,yerr=g_magerrs,color='seagreen',marker='.',ms=0.5,lw=0,elinewidth=0.5,label='green')
        
        axs[1].errorbar(c_r_dates,c_r_mags,yerr=c_r_magerrs,color='tomato',marker='.',ms=0.5,lw=0,elinewidth=0.5,label='red')
        axs[1].set_xlim(58430,58460)
        
        
        for i in range(2):
            axs[i].invert_yaxis()
            axs[i].set_xlabel('Date (MJD)')
            axs[i].set_ylabel('Mag')
            axs[i].legend(loc='lower right',fontsize='x-small')

        plt.subplots_adjust(hspace=0.5)
        plt.show()
        plt.savefig(out_temp.replace('+++',id),format='svg')
        plt.clf()
        plt.close()
        for i in range(2):
            axs[i].clear()
        bar.next()

    bar.finish()

ids = list(pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/key.csv')['ID'])

create_plots(IDs=['2MASS_J20521101+4408102'],file_temp='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++_*.tbl',out_temp='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/fixing_cluster_remover/plots/+++.svg')


    
    

    