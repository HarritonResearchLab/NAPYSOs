import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
from sklearn import linear_model
from os import path

coords = []
with open('/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/coordinates.txt','r') as f:
    for line in f:
        line = line.replace('\n','')
        coords.append(line)


def cmd_with_linear(coord,cmd_path,save_path):
    
    df = pd.read_csv(cmd_path)
    g_minus_rs = np.array(df['g-r'])
    g = np.array(df['g'])
    x_errs = np.array(df['g-r magerr'])
    y_errs = np.array(df['g magerr'])

    OLS = stats.linregress(g_minus_rs, g)

    TwoD_g_minus_rs = np.atleast_2d(g_minus_rs).T
    TwoD_g = np.atleast_2d(g).T
    ransac = linear_model.RANSACRegressor()
    ransac.fit(TwoD_g_minus_rs, TwoD_g)
    inlier_mask = ransac.inlier_mask_
    line_y_ransac = ransac.predict(TwoD_g_minus_rs)

    OLS_Inlier = stats.linregress(g_minus_rs[inlier_mask], g[inlier_mask])

    #Plot 
    sns.set_style('whitegrid')
    plt.title(coord+" CMD")
    plt.ylabel('g')
    plt.xlabel('g-r')
    plt.gca().invert_yaxis()
    plt.errorbar(g_minus_rs, g,xerr=x_errs,yerr=y_errs,lw=0,elinewidth=0.5,marker='.',color='cornflowerblue')
    
    OLS_slope = round(float(OLS[0]),3)
    OLS_rsq = round((OLS.rvalue * OLS.rvalue),3)
    plt.plot(g_minus_rs, OLS.intercept + OLS.slope * g_minus_rs, 'black',
             label=('OLS '+r'$r^2$: '+ str(OLS_rsq)+'; Slope: '+str(OLS_slope)))
    
    OLS_i_rsq = round((OLS_Inlier.rvalue * OLS_Inlier.rvalue),3)
    OLS_i_slope = round(float(ransac.estimator_.coef_),3)
    plt.plot(g_minus_rs[inlier_mask], OLS_Inlier.intercept + OLS_Inlier.slope * g_minus_rs[inlier_mask], 'red',
             label=('OLS Inlier '+r'$r^2$: '+ str(OLS_i_rsq)+'; Slope: '+str(OLS_i_slope)))
    
    plt.legend(loc='lower center', fontsize=6)

    file_path = save_path.replace('+++++',coord)
    plt.savefig(file_path) 
    plt.clf()

for coordinate in coords:
    csv_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/cmd_data/+++++.csv'.replace('+++++',coordinate)
    plot_path = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/for_email/cmd_plots/+++++.png'.replace('+++++',coordinate)
    cmd_with_linear(coord=coordinate,cmd_path=csv_file,save_path=plot_path)


"""
def cmd_with_linear(cmd_path):
    df = pd.read_csv(cmd_path)
    g_minus_rs = np.array(df['g-r'])
    g = np.array(df['g'])
    OLS = stats.linregress(g_minus_rs, g)

    TwoD_g_minus_rs = np.atleast_2d(g_minus_rs).T
    TwoD_g = np.atleast_2d(g).T
    ransac = linear_model.RANSACRegressor()
    ransac.fit(TwoD_g_minus_rs, TwoD_g)
    inlier_mask = ransac.inlier_mask_
    line_y_ransac = ransac.predict(TwoD_g_minus_rs)

    OLS_Inlier = stats.linregress(g_minus_rs[inlier_mask], g[inlier_mask])

    plt.figure()
    plt.title(path.basename(cmd_path) + " CMD")
    plt.ylabel('g')
    plt.xlabel('g-r')
    plt.gca().invert_yaxis()
    plt.scatter(g_minus_rs, g)
    plt.plot(g_minus_rs, OLS.intercept + OLS.slope * g_minus_rs, 'r',
             label=('OLS r^2: ' + str(OLS.rvalue * OLS.rvalue)))
    plt.plot(g_minus_rs, line_y_ransac, color='g', label='RANSAC')
    plt.plot(g_minus_rs[inlier_mask], OLS_Inlier.intercept + OLS_Inlier.slope * g_minus_rs[inlier_mask], 'y',
             label='OLS (inlier) r^2: ' + str(OLS_Inlier.rvalue * OLS_Inlier.rvalue))
    plt.legend(loc='lower center', borderpad=0.3, handlelength=0.7)

    return {"r_squared": OLS.rvalue * OLS.rvalue, "r_in_squared": OLS_Inlier.rvalue * OLS_Inlier.rvalue,
            "slope": OLS.slope, "in_slope": OLS_Inlier.slope}


sns.set_style('whitegrid')
cmd_with_linear('/Users/noah/PycharmProjects/milesGroup/data/cmd_data/20:51:15.14+44:18:17.4.csv')
plt.show()
"""