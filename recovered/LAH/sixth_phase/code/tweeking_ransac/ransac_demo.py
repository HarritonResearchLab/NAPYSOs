def ransac_demo(cmd_data): #Declaring a function here
    #Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from sklearn.linear_model import RANSACRegressor
    from scipy import stats

    #Action

    #Get data

    df = pd.read_csv(cmd_data) #read the data from the data file that you gave the function
    g_minus_rs = np.array(df['g-r'])
    g = np.array(df['g'])

    ransac = RANSACRegressor() #Define the ransac method

    TwoD_g_minus_rs = np.atleast_2d(g_minus_rs).T #See below
    TwoD_g = np.atleast_2d(g).T #See below
    
    # Note: so do you see how we did this np.atleast_2d().T on our x and y data arrays? 
    #       so I'm not entirely sure why we do this or what it means, but I know RANSAC needs that data in this
    #       format so I modified it as such. I don't think it's important to understand this (as I don't really tbh),
    #       the most important thing is to just get how RANSAC works overall. This step is necessary but not conceptually
    #       important.  

    ransac.fit(TwoD_g_minus_rs, TwoD_g) #Fit RANSAC on the 2D data 

    inlier_mask = ransac.inlier_mask_ #get array which tells you which values are inliers 

    outlier_mask = np.logical_not(inlier_mask) #based on inlier_mask, get array which tells you which values are outliers

    ransac_x = np.arange(g_minus_rs.min(),g_minus_rs.max())[:, np.newaxis] #Create a list of x values for ransac to make a line of off
    ransac_y = ransac.predict(ransac_x) #Based on the above x values, predict the ransac best fit line

    #Run ordinary least squares on RANSAC defined inliers to get slope and r^2 values: 

    ols = stats.linregress(g_minus_rs[inlier_mask], g[inlier_mask])

    ols_slope = ols[0] #get slope from linregress function
    ols_rsq = (ols[2])**2 #ols[2] returns r value, so squaring it with **2 returns r^2
    
    # Now that you have those values, if you modify this script to fit it on hundreds of CMDs, you could save those values
    # to a list or something for later analysis (I want to show you this today after you show Pietro how to do OLS on multiple files)


    #Plot
    plt.scatter(g_minus_rs[inlier_mask],g[inlier_mask],color='cornflowerblue',marker='.') #Create a scatter plot of the inlier values
    plt.scatter(g_minus_rs[outlier_mask],g[outlier_mask],color='indianred',marker='.') #Same, but for the outliers 
    plt.plot(ransac_x,ransac_y,color='black') #Plot the best fit RANSAC line
    plt.set_xlabel('g-r')
    plt.set_ylabel('g')
    plt.gca().invert_yaxis() #Invert the y-axis, because mags (yktv)
    plt.show()




#So now I'll show you how to use the function. I already have a data file I want to use defined below: 
data_file = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/cmd_data/2MASS_J20475581+4329019.csv'
#You'll have to change that data file to a file that actually exists on your computer. Then you give it to the function
# and tell the function to run like this: 

ransac_demo(cmd_data=data_file) #and the function will run, blue dots are inliers, red dots are outliers 