import numpy as np
from personalastropy.ysospy.handy_scripts import sortData
from personalastropy.ysospy.plotting_funcs import plotLightCurve
from astropy.io import ascii
from personalastropy.ysospy.handy_scripts import queryCoordSimbad

def plotting_demo():
    coordinates = []
    with open('/home/thaddaeus/FMU/HRL/LAH/third_phase/data/obsids.txt','r') as f:
        for line in f:
            line = line.replace('\n','')
            coordinates.append(line)

    for coordinate in coordinates:
        red_file = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_'+coordinate+'_r.tbl'
        green_file = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_'+coordinate+'_g.tbl'
        
        r = ascii.read(red_file,format='ipac',delimiter='|') 
        g = ascii.read(green_file,format='ipac',delimiter='|') 

        red_dates = np.array(r['hjd'])
        red_mags = np.array(r['mag'])
        srd = sortData(x=red_dates,y=red_mags)[0]
        srm = sortData(x=red_dates,y=red_mags)[1]

        green_dates = np.array(g['hjd'])
        green_mags = np.array(g['mag'])
        sgd = sortData(x=green_dates,y=green_mags)[0]
        sgm = sortData(x=green_dates,y=green_mags)[1]

        save_path = '/home/thaddaeus/FMU/HRL/LAH/third_phase/images/lightcurve_plots/++++++++++'.replace('++++++++++',coordinate)
        save_path = save_path+'.png'
        obsid = queryCoordSimbad(coordinate,5)
        plotLightCurve(x=[green_dates,red_dates],y=[green_mags,red_mags],colors=['green','red'],x_label='HJD',y_label='Mag',plot_title=(obsid+' Lightcurve'),line_labels=['Green Band','Red Band'],plot_type='scatter',out_type=save_path,error_arrays='na')
        
plotting_demo()
