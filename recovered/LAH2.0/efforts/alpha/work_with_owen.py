#Imports
from scipy.stats import skewnorm
import numpy as np
import pandas as pd

#Action
def return_nearest(array, value, corresponding_array):
    #Import(s)
    array = np.array(array)
    nearest_index = (np.abs(array - value)).argmin()
    return corresponding_array[nearest_index]

def collect_errors(ids,csv_temp):
    #Import(s)
    import pandas as pd
    import numpy as np

    #Action
    all_xs = []
    all_ys = []
    all_xerrs = []
    all_yerrs = []

    for id in ids:
        csv_file = csv_temp.replace('+++',id)
        df = pd.read_csv(csv_file)
        x = list(df['g-r'])
        y = list(df['g'])
        xerr = list(df['g-r magerr'])
        yerr = list(df['g magerr'])
        all_xs = all_xs + x
        all_ys = all_ys + y
        all_xerrs = all_xerrs + xerr
        all_yerrs = all_yerrs + yerr

    return np.array(all_xs),np.array(all_ys),np.array(all_xerrs),np.array(all_yerrs)


csv_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/cmd_data/+++.csv'
key_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/key.csv'
all_ids = list(pd.read_csv(key_path)['ID'])

cvs = collect_errors(ids=all_ids,csv_temp=csv_temp)

#ALL ERRORS 
all_x = cvs[0]
all_y = cvs[1]
all_xerr = cvs[2]
all_yerr = cvs[3]

#SIMULATE DATA

def simulate_data(slope,intercept):
    #Import(s)
    import matplotlib.pyplot as plt

    #Action
    x_vals = []
    while len(x_vals) < np.random.normal(45,10,1)+1:
        x_val = np.random.uniform(low=1,high=2)
        x_vals.append(x_val)
    
    xerrs = []
    for item in list(x_vals):
        xerrs.append(return_nearest(all_x,item,all_xerr))

    y_vals = (slope*np.array(x_vals)) + intercept
    yerrs = []
    
    for item in list(y_vals):
        yerrs.append(return_nearest(all_y,item,all_yerr))
    
    x_vals = list(np.random.normal(x_vals,(2*np.array(xerrs)),len(xerrs)))
    y_vals = list(np.random.normal(y_vals,yerrs,len(yerrs)))
    
    #Plot
    plt.rcParams['font.family'] = 'serif'
    plt.errorbar(x_vals,list(y_vals),xerr=xerrs,yerr=yerrs,lw=0,elinewidth=1,ms=1,color='cornflowerblue')
    plt.gca().invert_yaxis()
    plt.xlabel('g-r')
    plt.ylabel('g')
    plt.show()
    
simulate_data(slope=1.027,intercept=19)

