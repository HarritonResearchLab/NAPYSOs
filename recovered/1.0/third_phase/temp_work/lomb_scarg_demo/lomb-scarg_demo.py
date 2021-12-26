from personalastropy.ysospy.variability import lombScargle
from astropy.io import ascii
import numpy as np

coordinates = []
with open('/home/thaddaeus/FMU/HRL/LAH/third_phase/data/obsids.txt','r') as f:
    for line in f:
        line = line.replace('\n','')
        coordinates.append(line)

for coordinate in coordinates:
    red_file = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_'+coordinate+'_r.tbl'
    r = ascii.read(red_file,format='ipac',delimiter='|') 
    red_dates = np.array(r['hjd'])
    red_mags = np.array(r['mag'])
    if len(red_mags) >3:
        save_path = '/home/thaddaeus/FMU/HRL/LAH/third_phase/images/lomb_scarg_demo/++++++++++'.replace('++++++++++',coordinate)
        save_path = save_path+'.png'
        lombScargle(id=coordinate,x=red_dates,y=red_mags,min_period=2,max_period=400,out_type=save_path)