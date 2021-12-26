#Import(s)
from personalastropy.ysospy.interpolation import returnGoodRegions, cmdPrep
from scipy import interpolate
from personalastropy.ysospy.handy_scripts import sortData
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from scipy.stats.mstats import theilslopes
from sklearn import linear_model 
from astropy.io import ascii
import seaborn as sns

path_temp = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_20:53:15.62+43:44:22.8_=.tbl'
g_file = path_temp.replace('=','g')
r_file = path_temp.replace('=','r')

r=ascii.read(r_file,format='ipac',delimiter='|')
g=ascii.read(g_file,format='ipac',delimiter='|')

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
    

    #ols regression one
    X = np.array(interf[2]).reshape(-1,1)
    y = np.array(interf[1])

    lr = linear_model.LinearRegression()
    lr.fit(X,y)

    #Ransac regression
    ransac = linear_model.RANSACRegressor()
    ransac.fit(X,y)
    inlier_mask = ransac.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)
    
    #Predict data 
    line_X = np.arange(X.min(),X.max())[:, np.newaxis]
    line_y = lr.predict(line_X)
    line_y_ransac = ransac.predict(line_X)

    #Compare 
    print("est. coefficients (lin reg, RANSAC):")
    print(lr.coef_,ransac.estimator_.coef_)

    #plot

    plt.scatter(X[inlier_mask],y[inlier_mask],color='yellowgreen',marker='.',label='inliers')
    plt.scatter(X[outlier_mask],y[outlier_mask],color='gold',marker='.',label='outliers')
    plt.plot(line_X,line_y,color='navy',linewidth=2,label='OLS')
    plt.plot(line_X,line_y_ransac,color='cornflowerblue',linewidth=2,label='RANSAC')
    plt.gca().invert_yaxis()
    plt.xlabel('g-r')
    plt.ylabel('g')
    plt.title('Color-Mag Diagram')
    plt.legend(loc='lower right')
    plt.show()


