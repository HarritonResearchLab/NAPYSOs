import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import ascii
from matplotlib import rcParams
import seaborn as sns
from scipy import interpolate

#data files
green_path = './lc_203230.44+435741.9_g.tbl'#CHANGE THIS PATH!!
red_path = './lc_203230.44+435741.9_r.tbl' #CHANGE THIS PATH!!
#read the files into table objects
g = ascii.read(green_path,format='ipac',delimiter='|')
r = ascii.read(red_path,format='ipac',delimiter='|')

#read the relevant data columns into lists....idk if the np.array() part is 
#nessessary though lol

green_dates = list(np.array(g['hjd']))
green_mags = list(np.array(g['mag']))
red_dates = list(np.array(r['hjd']))
red_mags = list(np.array(r['mag']))

#order them by date...because I was having issues with polynomial fits
#when some of the data points were out of order so I ordered them here 
sgm = []
sgd = list(np.sort(np.array(green_dates)))

for elem in sgd:
    newIndex = sgd.index(elem)
    oldIndex = green_dates.index(elem)
    sgm.append(green_mags[oldIndex])

srm = []
srd = list(np.sort(np.array(red_dates)))

for elem in srd:
    newIndex = srd.index(elem)
    oldIndex = red_dates.index(elem)
    srm.append(red_mags[oldIndex])

#get dates that aren't in green band but are in red band to interpolate
new_sgd = []

for item in srd:
    if item not in sgd:
        if item < max(sgd) and item > min(sgd): #So it doesn't extrapolate outside the green domain
            new_sgd.append(item)

#turn all the lists into arrays
new_sgd = np.array(new_sgd)
sgm = np.array(sgm)
sgd = np.array(sgd)

srm =np.array(srm)
srd = np.array(srd)

#Make the x values a little less unwieldly 
sgd = sgd - (2.458*10**6)
srd = srd - (2.458*10**6)
new_sgd = new_sgd - (2.458*10**6)


#interpolate function
f = interpolate.interp1d(sgd,sgm,kind='slinear')
new_sgm = f(new_sgd)

#plot
sns.set_style('darkgrid')
rcParams['font.family']='Nimbus Roman'
plt.figure(figsize=(1,2))
plt.scatter(sgd,sgm,marker='o',s=2,color='forestgreen',label='Green Band')
plt.scatter(srd,srm,marker='o',s=2,color='red',label='Red Band')
plt.scatter(new_sgd,new_sgm,marker='x',color='blue',s=0.5,label='Interpolated Green Data')
plt.legend(loc='lower center',borderpad=0.3,handlelength=0.7,)
plt.xlabel('HJD — 2458000 [days]')
plt.ylabel('Mag')
plt.title('Real and Interpolated Data Demo')
plt.show()

