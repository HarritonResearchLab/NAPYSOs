import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import re
from astropy.coordinates import SkyCoord
from astropy import units as u

# Make plot
plt.rcParams['font.family']='serif'
plt.rcParams['figure.figsize']=(5, 5)
ax = plt.gca()

# harmonic lines
plt.axline(xy1=(0, 0), slope=1, color='grey', ls='--')
plt.axline(xy1=(0, 0), slope=0.5, color='grey', ls='--')
plt.axline(xy1=(0, 0), slope=2, color='grey', ls='--')

# Beat lines
beat_x = np.linspace(0.01, 25, 300)
beat_x1 = np.linspace(0.01, 0.94, 50) # due to asymtotes
beat_x2 = np.linspace(1.06, 15, 200)

# first
plt.plot(beat_x, 1/(1/beat_x+1), color='grey', lw=0.5)

# second
plt.plot(beat_x1, 1/(1/beat_x1-1), color='grey', lw=0.5)
plt.plot(beat_x2, 1/(1/beat_x2-1), color='grey', lw=0.5)

# third
plt.plot(beat_x1, 1/(beat_x1-1)+1, color='grey', lw=0.5)
plt.plot(beat_x2, 1/(beat_x2-1)+1, color='grey', lw=0.5)


## get periods and such 

our_results = pd.read_csv('./results/meta_data.csv')
our_ras = np.array(our_results['RA'])
our_decs = np.array(our_results['DEC'])
our_pers = np.array(our_results['period'])
our_primaries = np.array(our_results['primary_class'])

period_mask = np.where(our_primaries=='p')

# FROEBRICH FIRST!

# FIX: ADD COUNTERS!!!!
# FIX: add legend labels

froebrich_raw = './recovered/2.0/efforts/august21/aug23-29/froebrich/froebrich_raw_data.txt'

fro_periodic_counter = 0 
fro_timescale_counter = 0
with open(froebrich_raw, 'r') as f:
    for line in f: 
        line_list = re.sub(' +', ',', line).split(',')
        fro_ra = line_list[1] # used to be set to float
        fro_dec = line_list[2]
        fro_per = float(line_list[3])
        c2 = SkyCoord(fro_ra,fro_dec,frame='fk5',unit=(u.deg,u.deg)) # used to be converted to string here    
        for our_ra, our_dec, our_per, our_class in zip(our_ras, our_decs, our_pers, our_primaries):
            c1 = SkyCoord(str(our_ra),str(our_dec),frame='fk5',unit=(u.deg,u.deg))
            
            sep = c1.separation(c2).arcsecond
            
            if sep <= 2: 
                if our_class == 'p':
                    if fro_periodic_counter ==0: 
                        ax.scatter([our_per], [fro_per], color='lightsteelblue', edgecolors='black', linewidths=0.5, zorder=2.5, label='Froebrich et al.')
                        fro_periodic_counter+=1 
                    else: 
                        ax.scatter([our_per], [fro_per], color='lightsteelblue', edgecolors='black', linewidths=0.5, zorder=2.5)
                        fro_periodic_counter+=1 
                else: 
                    ax.scatter([our_per], [fro_per], color='lightsteelblue', edgecolors='black', linewidths=0.5, marker='s', zorder=2.5)
                    fro_timescale_counter+=1 

# BHARDWAJ SECOND
bhard_df = pd.read_csv('./recovered/2.0/efforts/kappa/Bhardwaj_comparison/Bhardwaj_results.csv')
bhard_ras = np.array(bhard_df['RA'])
bhard_decs = np.array(bhard_df['DEC'])
bhard_pers = np.array(bhard_df['PER'])

bhard_per_counter = 0
bhard_timescale_counter = 0


for bhard_ra, bhard_dec, bhard_per in zip(bhard_ras, bhard_decs, bhard_pers): 
    c2 = SkyCoord(str(bhard_ra),str(bhard_dec),frame='fk5',unit=(u.deg,u.deg)) # used to be converted to string here    
    for our_ra, our_dec, our_per, our_class in zip(our_ras, our_decs, our_pers, our_primaries):
        c1 = SkyCoord(str(our_ra),str(our_dec),frame='fk5',unit=(u.deg,u.deg))
        
        sep = c1.separation(c2).arcsecond
        
        if sep <= 2: 
            if our_class == 'p':
                if bhard_per_counter == 0: 
                    ax.scatter([our_per], [bhard_per], color='#408ee0', edgecolors='black', linewidths=0.5, zorder=2.5, label='Bhardwaj et al.')
                    bhard_per_counter +=1 
                else: 
                    ax.scatter([our_per], [bhard_per], color='#408ee0', edgecolors='black', linewidths=0.5, zorder=2.5)
                    bhard_per_counter +=1 
            else: 
                ax.scatter([our_per], [bhard_per], color='#408ee0', edgecolors='black', linewidths=0.5, marker='s', zorder=2.5) 
                bhard_timescale_counter+=1

# resume matching results
print('fro: ', fro_periodic_counter, fro_timescale_counter, fro_periodic_counter+fro_timescale_counter)
print('bhard: ', bhard_per_counter, bhard_timescale_counter, bhard_per_counter+bhard_timescale_counter)


# Resume with plot

plt.xlabel('Our Period (d)')
plt.ylabel('Literature Period (d)')

plt.xlim(left=0.0, right=14)
plt.ylim(bottom=0.0, top=14)
plt.tick_params(axis='both',which='both',direction='in')
plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
ax.xaxis.set_minor_locator(AutoMinorLocator())    
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.legend(shadow=False, edgecolor='black', fancybox=False)
plt.savefig('./referee_changes/redo_period_recovery/period_comparison.png', dpi=200, bbox_inches='tight')