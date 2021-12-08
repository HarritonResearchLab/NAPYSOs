import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
#speedOfSound1 = speed of sound when T = 15 K

speedOfLight = 29979200000 #cm/s
speedOfSound1 = 7763.7 #cm/s; when T = 15 Kelvin
maxTime = 3153600000 #number of seconds in 100 years

timeArray = np.linspace(0,maxTime,1000)
 
#plot
nimbusRoman = {'fontname':'Nimbus Roman'}

plt.style.use('seaborn-darkgrid')
plt.plot(timeArray,(timeArray*speedOfLight),linewidth=6,alpha=0.7,color='orange',label='Light Speed')
plt.plot(timeArray,(timeArray*speedOfSound1),linewidth=6,alpha=0.7,color='blue',label='Speed of Sound @ T = 15K')
plt.yscale('log')
plt.xscale('log')
plt.title('Figure 1.',**nimbusRoman)
plt.ylabel('Distance (cm)',**nimbusRoman)
plt.xlabel('Time (s)',**nimbusRoman)
plt.legend(prop={'family':'Nimbus Roman'})
plt.show()