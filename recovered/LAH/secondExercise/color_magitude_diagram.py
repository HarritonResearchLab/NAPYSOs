import matplotlib.pyplot as plt 
from astropy.io import ascii
from matplotlib import rcParams
import numpy as np
import seaborn as sns
from scipy import interpolate
from personalastropy.caltech.Caltech import sort_data, find_differences, calculate_peaks
import astropy.units as u
from astropy.timeseries import LombScargle
from scipy.interpolate import interp1d
from personalastropy.caltech.plot_data import Plot_data

#Get data
green_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_g.tbl'
red_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_r.tbl'
g = ascii.read(green_path,format='ipac',delimiter='|')
r = ascii.read(red_path,format='ipac',delimiter='|')

green_dates = np.array(g['hjd'])-(2.458*10**6)
green_mags = np.array(g['mag'])
red_dates = np.array(r['hjd'])-(2.458*10**6)
red_mags = np.array(r['mag'])

#sort the data
srd = sort_data(unsorteddates=red_dates,unsortedmags=red_mags)[0] #red dates
srm = sort_data(unsorteddates=red_dates,unsortedmags=red_mags)[1] #red mags
sgd = sort_data(unsorteddates=green_dates,unsortedmags=green_mags)[0] #green dates
sgm = sort_data(unsorteddates=green_dates,unsortedmags=green_mags)[1] #green mags

Plot_data(x=[srd,sgd],y=[srm,sgm],colors=['red','green'],x_label='HJD',y_label='Mag',plot_title='Test',line_labels=['Red','Green'],plot_type='scatter',out_type='show',error_arrays='N/A')

def remove_data(dates,mags):
    import numpy as np
    dateslist = list(dates)
    magslist = list(mags)
    for elem in dateslist:
        if 455.63 < elem < 455.72:
            elem_index = dateslist.index(elem)
            dateslist.remove(elem)
            magslist.remove(magslist[elem_index])
        elif 437.3 < elem < 437.9:
            elem_index = dateslist.index(elem)
            dateslist.remove(elem)
            magslist.remove(magslist[elem_index])
    return np.array(dateslist)
    return np.array(magslist)


def plot_xdistances():
    green_differences = find_differences(x=sgd)
    red_differences = find_differences(x=srd)

    sns.set_style('darkgrid')
    rcParams['font.family']='Nimbus Roman'
    
    fig, (ax1,ax2) = plt.subplots(1,2,sharey=True) #sharey ? 
    
    ax1.set_title('Green Band Observations')
    ax1.set_xlabel('Difference in observation date')
    ax1.set_ylabel('Frequency')
    mean_difference = round(np.mean(green_differences),4)
    differences_std = round(np.std(green_differences),3)
    histo_label = 'Mean: '+str(mean_difference)+'; SD: '+str(differences_std)
    ax1.hist(green_differences, log=True, label=histo_label,color='green',bins=30)
    ax1.legend()

    ax2.set_title('Red Band Observations')
    ax2.set_xlabel('Difference in observation date')
    mean_difference = round(np.mean(red_differences),4)
    differences_std = round(np.std(red_differences),3)
    histo_label = 'Mean: '+str(mean_difference)+'; SD: '+str(differences_std)
    ax2.hist(red_differences, log=True, label=histo_label,color='red',bins=30)
    ax2.legend()

    plt.show()

def plot_color_mag():

    #Interpolate for missing green data
    f = interpolate.interp1d(sgd,sgm,kind='slinear')
    new_sgd = []
    for item in srd:
        if item not in sgd:
            if item < max(sgd) and item > min(sgd):
                new_sgd.append(item)

    new_sgd=np.array(new_sgd)
    new_sgm = f(new_sgd)
    #match the sizes
    srm = srm[:-1] #size matching should be automated (later)

    #plot    
    sns.set_style('darkgrid')
    rcParams['font.family']='Nimbus Roman'
    plt.title('Color-Magnitude Diagram for ZTF18aaxykqu')
    plt.ylabel('G')
    plt.xlabel('(G-R)')
    plt.scatter((new_sgm-srm),new_sgm,c='black',s=2)
    plt.ylim(21.5,19)
    plt.xlim(2,4)
    plt.show()

    #best routine for this: 
    '''
    #Interpolate for missing green data
    f = interpolate.interp1d(sgd,sgm,kind='slinear')
    new_sgd = []
    for item in srd:
        if item not in sgd:
            if item < max(sgd) and item > min(sgd):
                new_sgd.append(item)

    new_sgd=np.array(new_sgd)
    new_sgm = f(new_sgd)
    #match the sizes
    srm = srm[:-1] #size matching should be automated (later)

    #plot    
    sns.set_style('darkgrid')
    rcParams['font.family']='Nimbus Roman'
    plt.title('Color-Magnitude Diagram for ZTF18aaxykqu')
    plt.ylabel('G')
    plt.xlabel('(G-R)')
    plt.scatter((new_sgm-srm),new_sgm,c='black',s=2)
    plt.ylim(21.5,19)
    plt.xlim(2,4)
    plt.show()
    '''
def analyze_variations(): 
    t_days = sgd*u.day 
    y_mags = sgm*u.mag

    
    frequency, power = LombScargle(t_days,y_mags).autopower(method='fastchi2')
    
    ls = LombScargle(t_days,y_mags)
    probabilities = [0.1, 0.05,0.01]
    fap_levels = ls.false_alarm_level(probabilities)
    
    #Plot
    plt.plot(frequency,power)
    fap_line_x = np.linspace(0.005,2.5,3)
    fap_line_y = np.full(3,fap_levels[2])
    plt.plot(fap_line_x,fap_line_y,c='grey',ls='--',label='0.01')
    plt.xlim(2.5,0.005)
    plt.xscale('log')
    plt.xlabel('Frequency (1/d)')
    plt.ylabel('Power')
    plt.legend(title='FAP')
    plt.title('Green Band Periodogram')
    plt.show()
    #print(g['magerr'])

def plot_interpolate(): #This is really messy 
        
    #Get data
    green_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_g.tbl'
    red_path = '/home/thaddaeus/FMU/HRL/LAH/secondExercise/ZTF18aaxykqu_light_curve_files/lc_203230.44+435741.9_r.tbl'
    g = ascii.read(green_path,format='ipac',delimiter='|')
    r = ascii.read(red_path,format='ipac',delimiter='|')

    green_dates = np.array(g['hjd'])-(2.458*10**6)
    green_mags = np.array(g['mag'])
    red_dates = np.array(r['hjd'])-(2.458*10**6)
    red_mags = np.array(r['mag'])

    #sort the data
    srd = sort_data(unsorteddates=red_dates,unsortedmags=red_mags)[0] #red dates
    srm = sort_data(unsorteddates=red_dates,unsortedmags=red_mags)[1] #red mags
    sgd = sort_data(unsorteddates=green_dates,unsortedmags=green_mags)[0] #green dates
    sgm = sort_data(unsorteddates=green_dates,unsortedmags=green_mags)[1] #green mags
    
    
    new_sgd = []
    sgm = list(sgm)
    srm = list(srm)
    sgd = list(sgd)
    srd = list(srd)

    for item in srd:
        if item not in sgd:
            if item < max(sgd) and item > min(sgd): #So it doesn't extrapolate outside the green domain
                new_sgd.append(item)
   
    #srm.remove(srm[688])

    sgm = np.array(sgm)
    sgd = np.array(sgd)

    srm =np.array(srm)
    srd = np.array(srd)

    f = interpolate.interp1d(sgd,sgm,kind='quadratic')
    new_sgm = f(new_sgd)

    #plot
    sns.set_style('darkgrid')
    rcParams['font.family']='Nimbus Roman'
    fig, (ax1) = plt.subplots(1)    
    ax1.scatter(sgd,sgm,marker='o',s=2,color='forestgreen',label='Green Band Data')
    ax1.scatter(new_sgd,new_sgm,marker='x',color='blue',s=2,alpha=0.4,label='Interpolated Data')
    ax1.scatter(srd,srm,marker='o',s=2,color='red',label='Red Band Data')
    ax1.set_xlabel('HJD — 2458000 [days]')
    ax1.set_ylabel('Mag')
    ax1.set_ylim(22,15.5)
    ax1.legend(borderpad=0.3,handlelength=0.7,)
    ax1.set_title('ZTF18aaxykqu Light Curves')
    plt.show()

def plot_quality(): 
    green_quals = g['catflags']
    red_quals = r['catflags']
    fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
    ax1.hist(green_quals,bins=6,range=(-1,5))
    ax1.set_xlabel('Green Data Quality')
    ax1.set_ylabel('Frequency')
    ax1.set_xlim(0,2)
    ax2.hist(red_quals,bins=6,range=(-1,5))
    ax2.set_xlabel('Red Data Quality')
    ax2.set_ylabel('Frequency')
    plt.show()


