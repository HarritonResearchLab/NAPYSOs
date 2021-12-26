def first_test():
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt

    def calc_redchi(fit,x,y,yerr,N,n_free):
        return 1.0/(N-n_free)*sum(((fit-y)/yerr)**2)

    x = np.linspace(-4,4,100)
    true_y = 0.5*np.sin(x)+17
    y = 0.5*np.sin(x)+np.random.normal(0,0.3,len(x))+17
    yerr = np.random.normal(0,0.4,len(x))

    residuals = true_y-y
    rss = sum(residuals**2)
    mean_val = np.mean(y)
    tss = sum((mean_val-y)**2)
    r_sq = 1-(rss/tss)
    print(r_sq)
    # Calculate red chi
    print(calc_redchi(fit=true_y,x=x,y=y,yerr=yerr,N=len(x),n_free=0))

    #plt.errorbar(x,y,yerr=yerr,lw=0,elinewidth=1)
    plt.scatter(x,y,s=2)
    plt.plot(x,true_y)
    plt.show()

first_test()