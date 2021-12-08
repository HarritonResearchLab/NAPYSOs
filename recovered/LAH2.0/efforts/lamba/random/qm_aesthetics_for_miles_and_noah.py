# Import(s)
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns

# Action

sns.color_palette()
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.gca().tick_params(axis='both',which='both',direction='in')

'''
# PLOTTING ORDER TO MATCH COLORS #
Okay so I'm not sure how to match the colors to the Q-M
plot and this plot for each class because I don't actually know
how to call the colors used in the palette. So if you guys can find 
their names that would be nice, or you can just plot them in this
order so they have the same colors. You can get the shorthand terms
for each class by checking the LaTeX version of the paper in section 
\sec 5.1. 
ORDER: 
1. Periodic
2. Quasi-P Sym
3. Stochastic 
4. Aperiodic Dipper
5. quasi periodic dipper
6. Burster
7. Long timescale
'''

plt.axes().yaxis.set_minor_locator(MultipleLocator(0.25))
plt.axes().xaxis.set_minor_locator(MultipleLocator(0.25))
plt.tick_params(axis='both',which='both',direction='in')
plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
plt.xlabel(r'$K_s-[22\mu$'+'m] (mag)')
plt.ylabel(r'$K_s-[12\mu$'+'m] (mag)')
plt.legend(loc='lower left',fontsize=6,fancybox=False,edgecolor='black',shadow=False,ncol=2)
plt.show()