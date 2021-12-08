#Import(s)
import pandas as pd 
import numpy as np
from astropy.io import ascii
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt

#Action

coord_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/coordinates.txt'
csv_temp = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/cmd_data/+++.csv'
g_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/lightcurves/lc_+++_g.tbl'
r_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/lightcurves/lc_+++_r.tbl'

with open(coord_file,'r') as f:
    for line in f:
        line = line.replace('\n','')
        
        #lightcurve related

        red_data = r_file.replace('+++',line)
        green_data = g_file.replace('+++',line)
        r = ascii.read(red_data,format='ipac',delimiter='|')
        g = ascii.read(green_data,format='ipac',delimiter='|')

        red_dates = r['mjd']
        red_mags = r['mag']
        red_magerrs = r['magerr']

        green_dates = g['mjd']
        green_mags = g['mag']
        green_magerrs = g['magerr']

        #cmd/regression related
        data_file = csv_temp.replace('+++',line)
        df = pd.read_csv(data_file)

        x_vals = np.array(df['g-r']) #g-r vals
        y_vals = np.array(df['g']) #g mags

        ols = stats.linregress(x_vals,y_vals)

        slope = ols[0]
        intercept = ols[1]
        r_val = ols[2]
        p_val = ols[3]

        #Create plot

        fig, axs = plt.subplots(2,1)
        sns.set_style('white')

        axs[0].errorbar(x=red_dates,y=red_mags,lw=0,elinewidth=0.5,marker='.',color='tomato',label='r band')
        axs[0].errorbar(x=green_dates,y=green_mags,lw=0,elinewidth=0.5,marker='.',color='seagreen',label='g band')
        axs[0].invert_yaxis()
        axs[0].legend(loc='lower right')
        axs[0].set_xlabel('Time (MJD)')
        axs[0].set_ylabel('Mag')

        axs[1].scatter(x=x_vals,y=y_vals,marker='.',color='cornflowerblue')
        
        min_x = min(x_vals)
        max_x = max(x_vals)

        ols_line_x = np.linspace(min_x,max_x,4)
        ols_line_y = (slope*ols_line_x)+intercept

        ols_label = r'$r^2$' + ': ' + str(r_val**2) + '\n' + 'slope: ' + str(slope)
        axs[1].plot(ols_line_x,ols_line_y,ls='--',color='black',label=ols_label)
        axs[1].legend(loc='upper right')


        plot_temp = '/home/thaddaeus/FMU/HRL/LAH/sixth_phase/plots/mjd_and_hjd_cmds/+++.png'
        plot_path = plot_temp.replace('+++',line)

        plt.savefig(plot_path)
        plt.clf()
        plt.close()

        #Save regression results
        
        r_sq = str(r_val**2)
        result_string = line + ',' + str(slope) + ',' + r_sq + '\n'

        results_file = '/home/thaddaeus/FMU/HRL/LAH/sixth_phase/plots/mjd_and_hjd_cmds/results.txt'

        with open(results_file,'a') as f2:
            f2.write(result_string)