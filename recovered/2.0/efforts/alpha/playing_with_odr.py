import time
def odr_demo(ground_truth,y_int,n,x_std,y_std):
    #Import(s)
    from scipy import odr
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import linregress

    
    #Action

    x = np.linspace(0,10,n)
    x_clean = np.linspace(0,10,n)
    xerr = np.abs(np.random.normal(0,x_std,n))
    x = np.random.normal(x,xerr,n)

    y = (ground_truth*x)+y_int

    yerr = np.abs(np.random.normal(0,y_std,n))
    y = np.random.normal(y,yerr)

    def line(x,a,b):
        y = (a*x)+b
        return y
    def odr_line(B,x):
        y = B[0]*x+B[1]
        return y
    def perform_odr(x,y,xerr,yerr):
        linear = odr.Model(odr_line)
        mydata = odr.Data(x,y,wd=1./xerr,we=1./yerr)
        myodr = odr.ODR(mydata,linear,beta0=[0,0])
        output = myodr.run()
        return output

    odr_regression = perform_odr(x=x,y=y,xerr=xerr,yerr=yerr)
   
    odr_slope = odr_regression.beta[0]
    odr_slope_error = odr_regression.sd_beta[0]
    
    error = 100*odr_slope_error/odr_slope #in percent

    ols = linregress(x,y)

    ols_slope = ols[0]
    ols_int = ols[1]
    ols_r_sq = ols[2]**2

    #Plot
    '''
    plt.errorbar(x,y,xerr,yerr,'o',alpha=0.7)
    plt.plot(x,(ground_truth*x+y_int),color='green',label=('Ground Truth slope: '+str(ground_truth))) #Ground truth
    plt.plot(x,line(x,odr_regression.beta[0],odr_regression.beta[1]),label=('ODR\nSlope: '+str(round(odr_regression.beta[0],3))+'\nError: '+str(round(error,3))+'%'),color='purple') #ODR
    plt.plot(x,(ols_slope*x+ols_int),color='red',label=('OLS slope:'+str(round(ols_slope,3))+'; r^2: '+str(round(ols_r_sq,3))))
    plt.legend(fontsize=6)

    plt.title('N: '+str(n)+'; x_std:'+ str(x_std)+'; y_std: '+str(y_std))
  
    plt.show()
    '''
    return ols_r_sq,error



r_sqs = []
errors = []

###ADJUST THESE BELOW!!!!! \/\/\/\/\/\/
ground_truth = 2
y_int = 3
n=10
x_std = 6
y_std = 6
num_sim = 10000

start_time = time.time()

for i in range(num_sim):
    odr_run = odr_demo(ground_truth=ground_truth,y_int=y_int,n=n,x_std=x_std,y_std=y_std)
    r_sqs.append(odr_run[0])
    errors.append(odr_run[1])

end_time = time.time()

delta = end_time-start_time

print('Fit time: '+str(round(delta,0))+' sec')


#Plot
import matplotlib.pyplot as plt
plt.scatter(r_sqs,errors,s=0.25)
plt.xlabel(r'$r^2$')
plt.ylabel('Error')
plt.axhline(y=10,color='red',ls='--')
#plt.axvline(x=.80,color='red',ls='--')
plt.ylim(0,100)
plt.title('n: '+str(n)+'; x_std: '+str(x_std)+'; y_std: '+str(y_std)+'; num simulations: '+str(num_sim)+'; Ground truth: '+str(ground_truth))
plt.show()