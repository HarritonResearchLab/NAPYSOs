#imports
from astropy.utils import data
import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import ascii
import seaborn as sns
import sklearn.gaussian_process as gp
from personalastropy.ysospy.handy_scripts import sortData

def gaussianModel(dataSetX, dataSetY):
    sns.set_style('darkgrid')
    dataX2d = np.atleast_2d(dataSetX).T
    dataXlin = np.atleast_2d(np.linspace(min(dataSetX), max(dataSetX), 10000)).T
    kernel = gp.kernels.ConstantKernel(1.0, (1e-1, 1e3)) * gp.kernels.RBF(10.0, (1e-3, 1e3))
    model = gp.GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=0.1, normalize_y=True)
    model.fit(dataX2d, dataSetY)
    return model

def gaussianPlot(dataSetX, dataSetY, gaussModel,Title, xLabel, yLabel):
    dataXlin = np.atleast_2d(np.linspace(min(dataSetX), max(dataSetX), 10000)).T
    mag_pred, sigma = gaussModel.predict(dataXlin, return_std=True)
    #plt.figure()
    plt.plot(dataSetX, dataSetY, 'r.', label='Observations')
    plt.plot(dataXlin, mag_pred, 'b-', label='Prediction')
    plt.fill(np.concatenate([dataXlin, dataXlin[::-1]]), np.concatenate([mag_pred - 1.9600 * sigma, (mag_pred + 1.9600 * sigma)[::-1]]), alpha=.5, fc='b', ec ='None', label='95% confidence interval')
    plt.gca().invert_yaxis()
    plt.show()

path_temp = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_+++++_=.tbl'
red_file = path_temp.replace('+++++','20:56:59.32+43:47:52.9')
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

func = sortData(x=red_dates,y=[red_mags,red_magerrs])
srd = func[0]
srm=func[1][0]
sre=func[1][1]
gfunc=sortData(x=green_dates,y=[green_mags,green_magerrs])
sgd = gfunc[0]
sgm=gfunc[1][0]
sge=gfunc[1][1]


sgdModel=gaussianModel(dataSetX=sgd,dataSetY=sgm)
gaussianPlot(dataSetX=sgd,dataSetY=sgm,gaussModel=sgdModel,Title='',xLabel='', yLabel='')