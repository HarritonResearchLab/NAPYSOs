import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
from sklearn import linear_model

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

    mad = stats.median_abs_deviation(g)

    OLS = stats.linregress(g_minus_rs, g)

    TwoD_g_minus_rs = np.atleast_2d(g_minus_rs).T
    TwoD_g = np.atleast_2d(g).T

    #resid_thresh = MAD
    ransac1 = linear_model.RANSACRegressor()
    ransac1.fit(TwoD_g_minus_rs, TwoD_g)
    inlier_mask1 = ransac1.inlier_mask_

    OLS_Inlier_1 = stats.linregress(g_minus_rs[inlier_mask1], g[inlier_mask1])
    
    #resid_thresh = 2 MAD
    ransac2 = linear_model.RANSACRegressor(residual_threshold=(mad*2))
    ransac2.fit(TwoD_g_minus_rs, TwoD_g)
    inlier_mask2 = ransac2.inlier_mask_

    OLS_Inlier_2 = stats.linregress(g_minus_rs[inlier_mask2], g[inlier_mask2])

    #resid_thresh = 0.5 MAD
    ransac50 = linear_model.RANSACRegressor(residual_threshold=(mad*0.5))
    ransac50.fit(TwoD_g_minus_rs, TwoD_g)
    inlier_mask50 = ransac50.inlier_mask_

    OLS_Inlier_50 = stats.linregress(g_minus_rs[inlier_mask50], g[inlier_mask50])

    #resid_thresh = 1.5 MAD
    ransac15 = linear_model.RANSACRegressor(residual_threshold=(mad*1.5))
    ransac15.fit(TwoD_g_minus_rs, TwoD_g)
    inlier_mask15 = ransac15.inlier_mask_

    OLS_Inlier_15 = stats.linregress(g_minus_rs[inlier_mask15], g[inlier_mask15])

    #Plot 
    sns.set_style('darkgrid')
    fig, axs = plt.subplots(2,2,sharex=True,sharey=True)
    OLS_slope = round(float(OLS[0]),3)
    OLS_rsq = round((OLS.rvalue * OLS.rvalue),3)
    
    #RT = MAD
    axs[0,0].errorbar(g_minus_rs, g,xerr=x_errs,yerr=y_errs,lw=0,elinewidth=0.5,marker='.',color='cornflowerblue')
    axs[0,0].plot(g_minus_rs, OLS.intercept + OLS.slope * g_minus_rs, 'black',
             label=('OLS '+r'$r^2$: '+ str(OLS_rsq)+'; Slope: '+str(OLS_slope)))

    OLS_i_rsq1 = round((OLS_Inlier_1.rvalue * OLS_Inlier_1.rvalue),3)
    OLS_i_slope1 = round(float(ransac1.estimator_.coef_),3)
    axs[0,0].plot(g_minus_rs[inlier_mask1], OLS_Inlier_1.intercept + OLS_Inlier_1.slope * g_minus_rs[inlier_mask1], 'red',
             label=('OLS Inlier R_T = MAD'+r'$r^2$: '+ str(OLS_i_rsq1)+'; Slope: '+str(OLS_i_slope1)))
    axs[0,0].legend(loc='lower center', fontsize=4)
    
    #RT = 1.5MAD
    axs[0,1].errorbar(g_minus_rs, g,xerr=x_errs,yerr=y_errs,lw=0,elinewidth=0.5,marker='.',color='cornflowerblue')
    axs[0,1].plot(g_minus_rs, OLS.intercept + OLS.slope * g_minus_rs, 'black',
             label=('OLS '+r'$r^2$: '+ str(OLS_rsq)+'; Slope: '+str(OLS_slope)))

    OLS_i_rsq15 = round((OLS_Inlier_15.rvalue * OLS_Inlier_15.rvalue),3)
    OLS_i_slope15 = round(float(ransac15.estimator_.coef_),3)
    axs[0,1].plot(g_minus_rs[inlier_mask15], OLS_Inlier_15.intercept + OLS_Inlier_15.slope * g_minus_rs[inlier_mask15], 'red',
             label=('OLS Inlier R_T = 1.5MAD'+r'$r^2$: '+ str(OLS_i_rsq15)+'; Slope: '+str(OLS_i_slope15)))
    axs[0,1].legend(loc='lower center', fontsize=4)

    #RT = 2MAD
    axs[1,0].errorbar(g_minus_rs, g,xerr=x_errs,yerr=y_errs,lw=0,elinewidth=0.5,marker='.',color='cornflowerblue')
    axs[1,0].plot(g_minus_rs, OLS.intercept + OLS.slope * g_minus_rs, 'black',
             label=('OLS '+r'$r^2$: '+ str(OLS_rsq)+'; Slope: '+str(OLS_slope)))

    OLS_i_rsq75 = round((OLS_Inlier_2.rvalue * OLS_Inlier_2.rvalue),3)
    OLS_i_slope75 = round(float(ransac2.estimator_.coef_),3)
    axs[1,0].plot(g_minus_rs[inlier_mask2], OLS_Inlier_2.intercept + OLS_Inlier_2.slope * g_minus_rs[inlier_mask2], 'red',
             label=('OLS Inlier R_T = 2MAD'+r'$r^2$: '+ str(OLS_i_rsq75)+'; Slope: '+str(OLS_i_slope75)))
    axs[1,0].legend(loc='lower center', fontsize=4)

    #RT = 0.5MAD
    axs[1,1].errorbar(g_minus_rs, g,xerr=x_errs,yerr=y_errs,lw=0,elinewidth=0.5,marker='.',color='cornflowerblue')
    axs[1,1].plot(g_minus_rs, OLS.intercept + OLS.slope * g_minus_rs, 'black',
             label=('OLS '+r'$r^2$: '+ str(OLS_rsq)+'; Slope: '+str(OLS_slope)))

    OLS_i_rsq50 = round((OLS_Inlier_50.rvalue * OLS_Inlier_50.rvalue),3)
    OLS_i_slope50 = round(float(ransac50.estimator_.coef_),3)
    axs[1,1].plot(g_minus_rs[inlier_mask50], OLS_Inlier_50.intercept + OLS_Inlier_50.slope * g_minus_rs[inlier_mask50], 'red',
             label=('OLS Inlier R_T = 0.5MAD'+r'$r^2$ '+ str(OLS_i_rsq50)+'; Slope: '+str(OLS_i_slope50)))
    axs[1,1].legend(loc='lower center', fontsize=4)

    #save plot
    file_path = save_path.replace('+++++',coord)
    plt.savefig(file_path) 
    plt.clf()
    plt.close()

def cmd_analysis(coord,cmd_path,save_path):
    
    #Action

    df = pd.read_csv(cmd_path)
    g_minus_rs = np.array(df['g-r'])
    g = np.array(df['g'])
    if len(g) > 20: 
    #x_errs = np.array(df['g-r magerr'])
    #y_errs = np.array(df['g magerr'])

        #Calculate some statistics on data
        mad_y = stats.median_abs_deviation(g)
        x_percentiles = np.percentile(g_minus_rs,[5,10,32,68,90,95])
        y_percentiles = np.percentile(g,[10,90])
        x_std = np.std(g_minus_rs,ddof=1)
        x_mean = np.mean(g_minus_rs)
        x_mad = stats.median_abs_deviation(g_minus_rs)
        y_std = np.std(g,ddof=1)
        y_mean = np.std(g)
        #change data

        #Ordinary OLS
        OLS = stats.linregress(g_minus_rs, g)
        OLS_slope = round(float(OLS[0]),3)
        OLS_rsq = round((OLS.rvalue * OLS.rvalue),3)

        TwoD_g_minus_rs = np.atleast_2d(g_minus_rs).T
        TwoD_g = np.atleast_2d(g).T

        #resid_thresh = MAD
        if x_percentiles[0]<(x_mean-x_mad) and x_percentiles[5] > (x_mean+x_mad):
            ransac = linear_model.RANSACRegressor(residual_threshold=2.5*x_std)
            plt.title('Condition met')
        else:
            ransac = linear_model.RANSACRegressor()
        ransac.fit(TwoD_g_minus_rs, TwoD_g)
        inlier_mask = ransac.inlier_mask_
        outlier_mask = np.logical_not(inlier_mask)
        ransac_x = np.arange(g_minus_rs.min(),g_minus_rs.max())[:, np.newaxis]
        ransac_y = ransac.predict(ransac_x)f

        #ransac statistics
        inlier_percent = 100*len(g[inlier_mask])/len(g)

        #OLS on inliers
        OLS_Inlier = stats.linregress(g_minus_rs[inlier_mask], g[inlier_mask])
        OLS_Inlier_rsq = round(((OLS_Inlier.rvalue)**2),3)
        OLS_Inlier_slope = round((OLS_Inlier.slope),3)
        OLS_Inlier_intercept = round((OLS_Inlier.intercept),3)

        #Plot 
        sns.set_style('darkgrid')
        plt.scatter(g_minus_rs[inlier_mask],g[inlier_mask],color='cornflowerblue',marker='.')
        plt.scatter(g_minus_rs[outlier_mask],g[outlier_mask],color='indianred',marker='.')
        #OLS_label = ()
        plt.plot(g_minus_rs,(OLS.intercept + OLS.slope * g_minus_rs))
        #RANSAC_label = () 
        plt.plot(ransac_x,ransac_y)
        
        #Inlier_OLS_label = ()
        plt.plot(g_minus_rs[inlier_mask],(OLS_Inlier_intercept+OLS_Inlier_slope*(g_minus_rs[inlier_mask])))

        #add x percentile lines
        for percentile in x_percentiles:
            plt.axvline(x=percentile,ls='--')
        for percentile in y_percentiles:
            plt.axhline(y=percentile,ls='--')
        
        for i in [-2,-1,1,2]:
            plt.axvline(x=(x_mean+(i*x_std)),ls='--',color='red')
        

        plt.gca().invert_yaxis()
        file_path = save_path.replace('+++++',coord)
        plt.savefig(file_path) 
        plt.clf()
        plt.close()
    




for coordinate in coords:
    csv_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/cmd_data/+++++.csv'.replace('+++++',coordinate)
    plot_path = '/home/thaddaeus/FMU/HRL/LAH/sixth_phase/plots/better_plots/+++++.png'.replace('+++++',coordinate)
    cmd_analysis(coord=coordinate,cmd_path=csv_file,save_path=plot_path)