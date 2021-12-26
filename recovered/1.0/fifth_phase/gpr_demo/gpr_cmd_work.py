#Import(s)
from astropy.utils import data
import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import ascii
import sklearn.gaussian_process as gp

def sortData(x,y):
    #Import(s)
    import numpy as np
    
    #Action
    
    unsorted_dates = list(x)
    
    sorted_dates = list(np.sort(unsorted_dates))
    sorted_arrays = []
    for item in y:
        unsorted_list = list(item)
        sorted_list = []
        for elem in sorted_dates:
            #newIndex = sorted_dates.index(elem)
            oldIndex = unsorted_dates.index(elem)
            sorted_list.append(unsorted_list[oldIndex])
        sorted_array = np.array(sorted_list)
        sorted_arrays.append(sorted_array)

    sorted_dates = np.array(sorted_dates)
    
    return sorted_dates, sorted_arrays

def gprModel(x, y):
    dataX2d = np.atleast_2d(x).T #x.T transposes x 
    kernel = gp.kernels.ConstantKernel(1.0, (1e-1, 1e3)) * gp.kernels.RBF(10.0, (1e-3, 1e3))
    model = gp.GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=0.1, normalize_y=True)
    model.fit(dataX2d, y)
    return model, model.kernel_.get_params()

def gprPlot(x, y, Model,Title, xLabel, yLabel):
    dataXlin = np.atleast_2d(np.linspace(min(x), max(x), 10000)).T
    mag_predictions, sigma = Model.predict(dataXlin, return_std=True)
    plt.plot(x, y, 'r.', label='Observations')
    plt.plot(dataXlin, mag_predictions, 'b-', label='Prediction')
    plt.fill(np.concatenate([dataXlin, dataXlin[::-1]]), np.concatenate([mag_predictions - 1.9600 * sigma, (mag_predictions + 1.9600 * sigma)[::-1]]), alpha=.5, fc='b', ec ='None', label='95% confidence interval')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.title(Title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.show()

g_file = './lc_20:53:15.62+43:44:22.8_g.tbl'

g=ascii.read(g_file,format='ipac',delimiter='|')
green_mags = g['mag']
green_dates=g['hjd']
green_magerrs=g['magerr']

# Below I sort the mags and magerr arrays to a sorted array of the dates because sometimes the arrays
# because I've seen that fitting routines (e.g. np.polyfit, etc.) get messed up when arrays aren't sorted 
# correctly.  
sort_func = sortData(x=green_dates,y=[green_mags,green_magerrs]) 
sgd = sort_func[0]
sgm=sort_func[1][0]
sge=sort_func[1][1]

g_gpr = gprModel(x=sgd,y=sgm)
print('Parameters:')
print(g_gpr[1])
gprPlot(x=sgd,y=sgm,Model=g_gpr[0],Title='',xLabel='Time (MJD)', yLabel='Mag')