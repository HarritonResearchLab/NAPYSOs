#Import(s)
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model 
from scipy.stats import linregress

#Action

#Load up data
df = pd.read_csv('./cmd_data/20:55:03.01+44:10:51.9.csv')
g_minus_rs = np.array(df['g-r'])
X = g_minus_rs.reshape(-1,1)
g_vals = np.array(df['g'])
y = g_vals
xerrs = np.array(df['g-r magerr'])
yerrs = np.array(df['g magerr'])


#Ransac regression
ransac = linear_model.RANSACRegressor()
ransac.fit(X,y)
inlier_mask = ransac.inlier_mask_
outlier_mask = np.logical_not(inlier_mask)

#Predict data 
line_X = np.arange(X.min(),X.max())[:, np.newaxis]
line_y_ransac = ransac.predict(line_X)

#Return values 
print('RANSAC SLOPE:')
print(ransac.estimator_.coef_)


print(line_X,line_y_ransac)
#Plot
plt.errorbar(X[inlier_mask],y[inlier_mask],xerr=xerrs[inlier_mask],yerr=yerrs[inlier_mask],color='cornflowerblue',elinewidth=0.5,linewidth=0,marker='o',ms=1,label='inliers')
plt.errorbar(X[outlier_mask],y[outlier_mask],xerr=xerrs[outlier_mask],yerr=yerrs[outlier_mask],color='palevioletred',elinewidth=0.5,linewidth=0,marker='o',ms=1,label='outliers')
plt.gca().invert_yaxis()
plt.xlabel('g-r')
plt.ylabel('g')
plt.title('20:55:03.01+44:10:51.9 CMD')
plt.legend(loc='lower right',fontsize=6)
plt.show()




