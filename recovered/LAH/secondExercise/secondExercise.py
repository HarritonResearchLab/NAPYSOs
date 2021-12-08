import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import ascii
from matplotlib import rcParams
import seaborn as sns
from scipy import interpolate
from scipy import stats
import scipy.stats as stats

##get Data
green_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_g.tbl'
red_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_r.tbl'
g = ascii.read(green_path,format='ipac',delimiter='|')
r = ascii.read(red_path,format='ipac',delimiter='|')

green_dates = list(np.array(g['hjd']))
green_mags = list(np.array(g['mag']))
red_dates = list(np.array(r['hjd']))
red_mags = list(np.array(r['mag']))

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

new_sgd = []


for item in srd:
    if item not in sgd:
        if item < max(sgd) and item > min(sgd): #So it doesn't extrapolate outside the green domain
            new_sgd.append(item)

new_sgd = np.array(new_sgd)
sgm = np.array(sgm)
sgd = np.array(sgd)

srm =np.array(srm)
srd = np.array(srd)

sgd = sgd - (2.458*10**6)
srd = srd - (2.458*10**6)
new_sgd = new_sgd - (2.458*10**6)

#smooth? 

from scipy.signal import savgol_filter
smoothed_red = savgol_filter(srm,3,2)

#interpolate function
f = interpolate.interp1d(sgd,sgm,kind='slinear')
new_sgm = f(new_sgd)

#find peaks
from scipy.signal import find_peaks
fakepeaks, _ = find_peaks(smoothed_red,width=5)

#find the 'real' peaks - "minimums" in scipy's eyes
negative_smoothed_red = smoothed_red*-1
realpeaks, _ = find_peaks(negative_smoothed_red,width=5)
peaks = np.concatenate((fakepeaks,realpeaks))
#plot

#global plot declarations
sns.set_style('darkgrid')
rcParams['font.family']='Nimbus Roman'
fig, (ax1) = plt.subplots(1) #change it to (ax1,ax2) and plt.subplots(2) to incl the histo
#fig.suptitle('Investigation')
fig.suptitle('ZTF18aaxykqu Light Curves')


#Plot observed and interpolated data
#ax1.scatter(sgd,sgm,marker='o',s=3,color='forestgreen',label='Green Band Data')
ax1.scatter(srd,srm,marker='o',s=3,color='red',label='Red Band Data')
#ax1.scatter(srd,smoothed_red,s=2,color='yellow',label='smoothed')
#ax1.scatter(new_sgd,new_sgm,marker='x',color='blue',s=0.5,label='Interpolated Green Data')
#ax1.legend(loc='lower center',borderpad=0.3,handlelength=0.7,) ###show this legend!
ax1.scatter(srd[peaks],smoothed_red[peaks],color='purple',s=5,marker='x')
ax1.set_xlabel('HJD — 2458000 [days]')
ax1.set_ylabel('Mag')
ax1.set_ylim(18.4,15.5) #22.5,15.5 for green too

#plot distribution of interpolated green band - red band mag values
srm = list(srm)
srm.remove(srm[688])
srm = np.array(srm)

#ax2.hist((new_sgm-srm))

plt.show()



