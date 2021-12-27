# Import(s)
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator


# Action

df = pd.read_csv(r'C:\Users\Research\Documents\GitHub\NAPYSOs\recovered\2.0\efforts\september21\SEPTEMBER_FINAL_RESULTS\merged_results.csv')
periods = np.array(df['PER'])
class_mask = np.where(np.array(df['primary_class'])=='p')
periods = periods[class_mask]
timescales = np.array(df['PER'])
# Plot
plt.rcParams['font.family']='serif'
fig, ax = plt.subplots()
ax.set_box_aspect(1)
#ax.set_aspect('equal')

ax.hist(timescales, color='lightsteelblue', edgecolor='black')

#plt.hist(periods,range=[0,20],bins=16,color='#0051a2',edgecolor='black')
ax.set(xlabel='Period (d)', ylabel='Frequency')
    
plt.tick_params(axis='both',which='both',direction='in')
plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
ax.xaxis.set_minor_locator(MultipleLocator(25))
plt.yscale('log')
#plt.minorticks_on()
plt.show()
#plt.savefig('timescales.png', dpi=250, bbox_inches='tight')