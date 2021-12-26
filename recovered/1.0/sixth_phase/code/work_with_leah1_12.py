#Import(s)
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii 
import pandas as pd

#Action

coord_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/coordinates.txt'
g_temp = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/lightcurves/lc_+++_g.tbl'
r_temp = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/lightcurves/lc_+++_r.tbl'
csv_temp = ''

with open(coord_file,'r') as f:
    for line in f:
        line = line.replace('\n','')

        red_data = r_temp.replace('+++',line)
        green_data = g_temp.replace('+++',line)

        r = ascii.read(red_data,format='ipac',delimiter='|')
        g = ascii.read(green_data,format='ipac',delimiter='|')

        red_dates = r['mjd']
        red_mags = r['mag']

        green_dates = g['mjd']
        green_mags = g['mag']

        #Create a plot

        fig, axs = plt.subplots(2,1)

        axs[0].scatter(x=red_dates,y=red_mags,s=3,marker='o',color='tomato',label='Red')
        axs[0].scatter(x=green_dates,y=green_mags,s=3,marker='o',color='seagreen',label='Green')
        axs[0].set_ylabel('Mag')
        axs[0].set_xlabel('Time (MJD)')
        axs[0].invert_yaxis()
        axs[0].legend(loc='lower right')

        

        plt.subplots_adjust(hspace=0.35)
        plot_path = '/home/thaddaeus/lc_plots/+++.png'.replace('+++',line)
        plt.savefig(plot_path)
        plt.clf()
        plt.close()




