#Import(s)
from astropy.io import ascii
from personalastropy.ysospy.interpolation import cmdPrep2, returnGoodRegions
from personalastropy.ysospy.handy_scripts import sortData
import pandas as pd
import numpy as np

#Action
coords = []
with open('coordinates.txt','r') as f:
    for line in f:
        line = line.replace('\n','')
        coords.append(line)


for coord in coords:
    g_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/lightcurves/lc_+++++_g.tbl'.replace('+++++',coord)
    r_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/lightcurves/lc_+++++_r.tbl'.replace('+++++',coord)

    g = ascii.read(g_file,format='ipac',delimiter='|')
    r = ascii.read(r_file,format='ipac',delimiter='|')

    green_mags = g['mag']
    green_dates=g['hjd']
    green_magerrs=g['magerr']

    red_mags=r['mag']
    red_dates=r['hjd']
    red_magerrs=r['magerr']

    if len(red_mags)>20 and len(green_mags)>20:
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


        cmdData = cmdPrep2(x=sgd,y=sgm,magerrs=sge,gdis=green_date_intervals,x2=srd,y2=srm,magerrs2=sre,kind='linear')
        
        csv_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/cmd_data/+++++.csv'.replace('+++++',coord)
        col_1 = list(cmdData[2])
        col_2 = list(cmdData[1])
        col_3 = list(cmdData[3])
        col_4 = list(cmdData[4])
        
        zipped_list = list(zip(col_1,col_2,col_3,col_4))

        df = pd.DataFrame(zipped_list,columns=['g-r','g','g-r magerr','g magerr'])    
        
        #df.to_csv(csv_file,index=False)
    else:
        print(coord)




