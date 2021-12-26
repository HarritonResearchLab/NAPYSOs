import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

df = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv')

nus = np.array(df['NU'])
timescales = np.array(df['PER'])
classes = np.array(df['primary_class'])

 # Split into classes
def return_arrays(class_name):
    arrays_arr = [nus, timescales]
    return [i[np.where(classes==class_name)] for i in arrays_arr]

p_nus, p_timescales = return_arrays('p')
mp_nus, mp_timescales = return_arrays('mp')
qpd_nus, qpd_timescales = return_arrays('qpd')
apd_nus, apd_timescales = return_arrays('apd')
qps_nus, qps_timescales = return_arrays('qps')
b_nus, b_timescales = return_arrays('b')
l_nus, l_timescales = return_arrays('l')
s_nus, s_timescales = return_arrays('s')

# Initial settings
plt.rcParams['font.family'] = 'serif'
plt.gca().tick_params(axis='both',which='both',direction='in')

# Plot classes
plt.scatter(p_nus,p_timescales,marker='o',color='C0',label='Periodic',edgecolor='black',linewidths=0.5)

# Quasi sym 
plt.scatter(qps_nus,qps_timescales,marker='o', color='C1', label='Quasi-periodic',edgecolor='black',linewidths=0.5)

# Stochastic
plt.scatter(s_nus,s_timescales,marker='o', color='C2',label='Stochastic',edgecolor='black',linewidths=0.5)

# Aper dipper
plt.scatter(apd_nus,apd_timescales,marker='o', color='C3', label='Aperiodic dipper',edgecolor='black',linewidths=0.5)

# Quasi periodic dippers
plt.scatter(qpd_nus,qpd_timescales,marker='o', color='C4', label='Quasi-periodic dipper',edgecolor='black',linewidths=0.5)

# Bursters
plt.scatter(b_nus,b_timescales,marker='o', color='C5', label='Burster',edgecolor='black',linewidths=0.5)

# Long timescales
plt.scatter(l_nus,l_timescales,marker='o', color='black', label='Long-timescale', edgecolor='black',linewidths=0.5)

# Multi periodics
plt.scatter(mp_nus, mp_timescales, marker='o', color='C6', label='Multi-periodic', edgecolor='black',linewidths=0.5)

# Final settings
plt.xlabel('Nu')
plt.ylabel('Timescale (d)')
plt.yscale('log')
#plt.axes().yaxis.set_minor_locator(MultipleLocator(0.05))
plt.axes().xaxis.set_minor_locator(AutoMinorLocator())
plt.tick_params(axis='both',which='both',direction='in')
plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)

lgnd = plt.legend(loc='center right',fontsize=6,fancybox=False,edgecolor='black',shadow=False,ncol=2)

for i in range(8): 
    lgnd.legendHandles[i]._sizes = [18]

plt.show()
