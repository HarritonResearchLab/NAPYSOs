import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from dust_extinction.parameter_averages import CCM89, F99
from matplotlib.ticker import AutoMinorLocator

# Create wavelengths array.
wav = np.arange(0.4, 0.77, 0.01)*u.micron

g_wav = 0.4723*u.micron
r_wav = 0.634*u.micron

# Calculate and Plot
#plt.style.use('seaborn-darkgrid')
plt.rcParams['font.family']='serif'
plt.rcParams['mathtext.fontset']='dejavuserif'
fig, axs = plt.subplots(2,2)

# Clayton, Cardelli, & Mathis (1989)
for model in [CCM89]:#, F99]:
    for R in (2.0,3.0,4.0):
        # Initialize the extinction model
        ext = model(Rv=R)
        axs[0,0].plot(1/wav, ext(wav), label='R='+str(R))
        axs[0,0].scatter(1/g_wav,ext(g_wav),color='green',marker='+')
        axs[0,0].scatter(1/r_wav,ext(r_wav),color='red',marker='+')

axs[0,0].set_xlabel('$\lambda^{-1}$ ($\mu$m$^{-1}$)')
axs[0,0].set_ylabel('A($\lambda$) / A(V)')
axs[0,0].legend(loc='best',fancybox=False,edgecolor='black',shadow=False,fontsize=7)
axs[0,0].tick_params(axis='both',which='both',direction='in')
axs[0,0].tick_params(which='both',bottom=True,top=True,left=True,right=True)
axs[0,0].tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].set_title('CCM89 Extinction Law')

# Fitzpatrick (1999)
for model in [F99]:
    for R in (2.0,3.0,4.0):
        # Initialize the extinction model
        ext = model(Rv=R)
        axs[0,1].plot(1/wav, ext(wav), label='R='+str(R))
        axs[0,1].scatter(1/g_wav,ext(g_wav),color='green',marker='+')
        axs[0,1].scatter(1/r_wav,ext(r_wav),color='red',marker='+')

axs[0,1].set_xlabel('$\lambda^{-1}$ ($\mu$m$^{-1}$)')
axs[0,1].set_ylabel('A($\lambda$) / A(V)')
axs[0,1].legend(loc='best',fancybox=False,edgecolor='black',shadow=False,fontsize=7)
axs[0,1].tick_params(axis='both',which='both',direction='in')
axs[0,1].tick_params(which='both',bottom=True,top=True,left=True,right=True)
axs[0,1].tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].set_title('F99 Extinction Law')

axs[1,0].axis('off')
axs[1,1].axis('off')

plt.subplots_adjust(wspace=0.5)

plt.show()