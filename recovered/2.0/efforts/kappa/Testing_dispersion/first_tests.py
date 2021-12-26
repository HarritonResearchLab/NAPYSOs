def first_test(id,pow,q,fap,data_temp,period,plot_path): 
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator

    # Action

    df = pd.read_csv(data_temp.replace('+++',id))
    mjds = np.array(df['mjd'])
    mags = np.array(df['mag'])
    magerrs = np.array(df['magerr'])

    # Get stds from raw light curve
    raw_stds = []
    for i in range(0,81,20):
        lp = i
        up = i+20
        date_mask = np.logical_and(mjds<np.percentile(mjds,up),mjds>=np.percentile(mjds,lp))
        masked_mags = mags[date_mask]
        std = np.std(masked_mags)
        raw_stds.append(float(std))

    # Fold the light curve
    phased_dates = (np.mod(mjds, period))/period
    
    # Get stds from folded light curve
    folded_stds = []
    for i in range(0,81,20):
        lp = i
        up = i+20
        date_mask = np.logical_and(phased_dates<np.percentile(phased_dates,up),phased_dates>=np.percentile(phased_dates,lp))
        masked_mags = mags[date_mask]
        std = np.std(masked_mags)
        folded_stds.append(float(std))

    raw_stds = np.array(raw_stds)
    folded_stds = np.array(folded_stds)
    std_division = folded_stds/raw_stds
    mean_std_division = np.mean(std_division)

    raw_variances = raw_stds**2
    folded_variances = folded_stds**2
    variance_division = folded_variances/raw_variances
    mean_variance_division = np.mean(variance_division)
    # Plot
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'


    fig, axs = plt.subplots(2,2)

    axs[0,0].set_title('FLC; Q: '+str(round(q,3)),fontsize='medium')
    axs[0,0].scatter(phased_dates,mags,color='#408ee0',marker='o',s=2)
    axs[0,0].scatter(phased_dates+1,mags,color='#408ee0',marker='o',s=2)
    xlabel = 'Phase (P = '+str(round(period, 3))+' d)'
    axs[0,0].set_xlabel(xlabel)
    axs[0,0].set_ylabel('Mag')
    axs[0,0].invert_yaxis()
    axs[0,0].locator_params(axis='y', nbins=5)
    axs[0,0].tick_params(axis='both',which='both',direction='in')
    axs[0,0].tick_params(which='both',bottom=True,top=True,left=True,right=True)
    axs[0,0].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
    axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())

    axs[0,1].set_title('Mean quotient: '+str(round(mean_std_division,3)),fontsize='medium')
    xs = np.array(range(1,6))
    axs[0,1].plot(xs,raw_stds,label='Light curve',color='#408ee0',ls='solid')
    axs[0,1].plot(xs,folded_stds,label='Folded light curve',color='#89bff8',ls='dashed')
    axs[0,1].scatter(xs,raw_stds,color='#408ee0',marker='o')
    axs[0,1].scatter(xs,folded_stds,color='#89bff8',marker='s')
    axs[0,1].set_ylabel(r'$\sigma$')
    axs[0,1].set_xlabel('Segment')
    axs[0,1].tick_params(axis='both',which='both',direction='in')
    axs[0,1].tick_params(which='both',bottom=True,top=True,left=True,right=True)
    axs[0,1].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    axs[0,1].set_xticks(range(1,len(xs)+1))
    axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
    axs[0,1].legend(loc='lower right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)
    
    axs[1,1].set_title('Mean quotient: '+str(round(mean_variance_division,3)),fontsize='medium')
    xs = np.array(range(1,6))
    axs[1,1].plot(xs,raw_variances,label='Light curve',color='#408ee0',ls='solid')
    axs[1,1].plot(xs,folded_variances,label='Folded light curve',color='#89bff8',ls='dashed')
    axs[1,1].scatter(xs,raw_variances,color='#408ee0',marker='o')
    axs[1,1].scatter(xs,folded_variances,color='#89bff8',marker='s')
    axs[1,1].set_ylabel(r'$\sigma^2$')
    axs[1,1].set_xlabel('Segment')
    axs[1,1].tick_params(axis='both',which='both',direction='in')
    axs[1,1].tick_params(which='both',bottom=True,top=True,left=True,right=True)
    axs[1,1].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    axs[1,1].set_xticks(range(1,len(xs)+1))
    axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())
    axs[1,1].legend(loc='lower right',fontsize=5,fancybox=False,edgecolor='black',shadow=False)

    axs[1,0].errorbar(mjds, mags, yerr=magerrs, lw=0,elinewidth=0.5)
    axs[1,0].scatter(mjds,mags,s=2)
    axs[1,0].invert_yaxis()
    axs[1,0].set_ylabel('Mag')
    axs[1,0].set_xlabel('MJD')
    axs[1,0].set_title('Light Curve', fontsize=10)
    axs[1,0].locator_params(axis='y', nbins=5)
    axs[1,0].tick_params(axis='both',which='both',direction='in')
    axs[1,0].tick_params(which='both',bottom=True,top=True,left=True,right=True)
    axs[1,0].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
    axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
    axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())

    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    plt.savefig(plot_path.replace('+++',id))
    plt.clf()
    plt.close()
    
    # Return 
    return id+','+str(pow)+','+str(period)+','+str(fap)+','+str(mean_variance_division)+','+str(q)

### RUN IT
import pandas as pd
key_df = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results.csv')
ids = list(key_df['ID'])
periods = list(key_df['PER1'])
faps = list(key_df['99%_FAP'])
powers = list(key_df['POW1'])
qs = list(key_df['Q'])

data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
plot_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Testing_dispersion/first_plots/+++.png'

with open('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Testing_dispersion/results.csv','a') as f:
    f.write('ID,POW,PER,FAP,QUOTIENT,Q'+'\n')
    for id, per, pow, fap, q in zip(ids,periods, powers, faps, qs): 
        f.write(first_test(id=id,pow=pow,q=q,fap=fap, data_temp=data_temp,period=per,plot_path=plot_path)+'\n')

    


    
