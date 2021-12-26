# Import(s)
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Action

df = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv')
periods = np.array(df['PER'])
class_mask = np.where(np.array(df['primary_class'])=='p')
periods = periods[class_mask]

# Plot
plt.rcParams['font.family']='serif'
fig, ax = plt.subplots()
ax.set_box_aspect(1)
#ax.set_aspect('equal')

ax.hist(periods, color='#408ee0',edgecolor='black')
#plt.hist(periods,range=[0,20],bins=16,color='#0051a2',edgecolor='black')
ax.set(xlabel='True Rotation Period (d)', ylabel='Frequency')
    
plt.tick_params(axis='both',which='both',direction='in')
plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
plt.minorticks_on()
plt.savefig('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/september21/27-30/for lah/rotation_periods.png', dpi=250, bbox_inches='tight')