# Import(s)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Action
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
fig, axs = plt.subplots(2,2)
plt.subplots_adjust(hspace=0.4,wspace=0.5)

# Plot

# axs[0,0]

x1 = np.random.normal(0,1,1000)
axs[0,0].hist(x1,color='#408ee0',label='Histogram')

axs[0,0].set_xlabel('x')
axs[0,0].set_ylabel('Frequency')

axs[0,0].legend(loc='upper left',fontsize=6,fancybox=False,edgecolor='black',shadow=False)
axs[0,0].tick_params(axis='both',which='both',direction='in')
axs[0,0].tick_params(bottom=True,top=True,left=True,right=True)
axs[0,0].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())

# axs[0,1]
x2 = np.linspace(0,5,5)
y1 = x2
y2 = 2*x2
y3 = 3*x2

axs[0,1].scatter(x2,y1,color='#408ee0',marker='o',label=r'$\frac{\Delta y}{\Delta x}$ = 1')
axs[0,1].scatter(x2,y2,color='#89bff8',marker='s',label=r'$\frac{\Delta y}{\Delta x}$ = 2')
axs[0,1].scatter(x2,y3,color='#0051a2',marker='d',label=r'$\frac{\Delta y}{\Delta x}$ = 3')

axs[0,1].set_xlabel('x')
axs[0,1].set_ylabel('y')

axs[0,1].legend(loc='upper left',fontsize=7,fancybox=False,edgecolor='black',shadow=False)
axs[0,1].tick_params(axis='both',which='both',direction='in')
axs[0,1].tick_params(bottom=True,top=True,left=True,right=True)
axs[0,1].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())

# axs[1,0]
axs[1,0].plot(x2,y1,color='#408ee0',ls='solid',label=r'$\frac{\Delta y}{\Delta x}$ = 1')
axs[1,0].plot(x2,y2,color='#89bff8',ls='dashed',label=r'$\frac{\Delta y}{\Delta x}$ = 2')
axs[1,0].plot(x2,y3,color='#0051a2',ls='dotted',label=r'$\frac{\Delta y}{\Delta x}$ = 3')

axs[1,0].set_xlabel('x')
axs[1,0].set_ylabel('y')

axs[1,0].legend(loc='upper left',fontsize=7,fancybox=False,edgecolor='black',shadow=False)
axs[1,0].tick_params(axis='both',which='both',direction='in')
axs[1,0].tick_params(bottom=True,top=True,left=True,right=True)
axs[1,0].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())

# axs[1,1]
x3 = np.linspace(-np.pi,np.pi,30)
y3 = np.sin(x3)

axs[1,1].scatter(x3,y3,color='#408ee0',marker='o',label=r'$\sin(x)$')

axs[1,1].set_xlabel('x')
axs[1,1].set_ylabel('y')

axs[1,1].legend(loc='upper left',fontsize=7,fancybox=False,edgecolor='black',shadow=False)
axs[1,1].tick_params(axis='both',which='both',direction='in')
axs[1,1].tick_params(bottom=True,top=True,left=True,right=True)
axs[1,1].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
axs[1,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())

# Show plot
plt.show()