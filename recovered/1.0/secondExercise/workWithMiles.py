from astropy.utils import data
import matplotlib.pyplot as plt 
from astropy.io import ascii
from matplotlib import rcParams
import numpy as np
import seaborn as sns
from scipy.signal import find_peaks

#get data
green_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_g.tbl'
red_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_r.tbl'
g = ascii.read(green_path,format='ipac',delimiter='|')
r = ascii.read(red_path,format='ipac',delimiter='|')

green_dates = np.array(g['hjd'])-(2.458*10**6)
green_mags = np.array(g['mag'])
red_dates = np.array(r['hjd'])-(2.458*10**6)
red_mags = np.array(r['mag'])


def sort_data(unsorteddates,unsortedmags):
    #Import
    import numpy as np
    #Action
    unsorteddates= list(unsorteddates)
    unsortedmags=list(unsortedmags)
    
    sorteddates=list(np.sort(unsorteddates))
    sortedmags = []
    for elem in sorteddates:
        newIndex = sorteddates.index(elem)
        oldIndex = unsorteddates.index(elem)
        sortedmags.append(unsortedmags[oldIndex])
    
    sorteddates = np.array(sorteddates)
    sortedmags = np.array(sortedmags)
    
    return sorteddates, sortedmags 

srd = sort_data(unsorteddates=red_dates,unsortedmags=red_mags)[0]
srm = sort_data(unsorteddates=red_dates,unsortedmags=red_mags)[1]

#smooth data
from scipy.signal import savgol_filter
#smoothed_red = savgol_filter(srm,3,2)

#Find peaks function
def calculate_peaks(data_array,widthval):
    fakepeaks, _ = find_peaks(data_array,width=widthval)
    negative_array = data_array*-1
    realpeaks, _ = find_peaks(negative_array,width=widthval)
    outarray = np.concatenate((fakepeaks,realpeaks))
    return(outarray)

peaks = calculate_peaks(data_array=srm,widthval=4)

#plot
sns.set_style('darkgrid')
rcParams['font.family']='Nimbus Roman'
plt.title('Investigating ZTF18aaxykqu')
plt.ylabel('Mag')
plt.xlabel('HJD — 2458000 [days]')
plt.scatter(red_dates,red_mags,c='red',s=2,label='Red Band Data')
plt.scatter(srd[peaks],srm[peaks],c='purple',s=2,label='Relative extrema')
plt.legend()
plt.ylim(18.5,15.5)
plt.show()
