# Import(s)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

first_df = pd.read_csv('/home/thaddaeus/Keppler Q.csv')
second_df = pd.read_csv('/home/thaddaeus/Downsampled Keppler Q.csv')

for x_column in first_df.columns.values.tolist(): 
    for y_column in second_df.columns.values.tolist():
        plt.scatter(first_df[x_column], second_df[y_column])
        plt.xlabel(x_column)
        plt.ylabel(y_column)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.savefig('/home/thaddaeus/plots/'+x_column+':'+y_column+'.png') 
        plt.clf()
        plt.close()