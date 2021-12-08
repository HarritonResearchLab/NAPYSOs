#Import(s)
import pandas as pd
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

#Action

plt.title(r'$\pi^x-\int_{0}^{2}x^2dx+\frac{x^2}{\sqrt{3x-2}}$')
plt.show()

'''
MATPLOTLIB: Plots 
    Plotting functions
        plt.plot(x,y) #line
        plt.scatter(x,y) #scatter plot
        plt.errorbar(x,y,xerr,yerr)
        plt.hist(x)

    #format functions
    plt.xlabel() #set xlabel, e.g. plt.xlabel('Time')
    plt.ylabel() #set ylabel, e.g. 
    plt.title() #change title, e.g. plt.title('Cool plot')
    plt.legend()
    plt.gca().invert_yaxis() #invert yaxis
    plt.show() #Show the plot
    plt.savefig('path') #Save plot to location in string
    plt.close() #close plot if you're making a lot of plots

'''
