def match_to_file(our_results, froebrich_raw, outfile):
    # Import(s)
    import numpy as np
    import pandas as pd
    import re
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    # Action
    fro_ras = np.array([])
    fro_decs = np.array([])
    fro_pers = np.array([])
    
    with open(froebrich_raw, 'r') as f:
        for line in f: 
            line_list = re.sub(' +', ',', line).split(',')
            fro_ras = np.append(fro_ras, float(line_list[1]))
            fro_decs = np.append(fro_decs, float(line_list[2]))
            fro_pers = np.append(fro_pers, float(line_list[3]))
            
    results_df = pd.read_csv(our_results)
    classes = np.array(results_df['primary_class'])
    
    our_ids = np.array(results_df['ID'])
    our_ras = np.array(results_df['RA'])
    our_decs = np.array(results_df['DEC'])
    our_pers = np.array(results_df['period'])
    
    shared_ids = np.array([])
    shared_our_pers = np.array([])
    shared_fro_pers = np.array([])
    shared_classes = np.array([])

    counter = 0 
    
    for fro_ra, fro_dec, fro_per in zip(fro_ras, 
                                        fro_decs, 
                                        fro_pers):
        num = 0
        
        c2 = SkyCoord(str(fro_ra),str(fro_dec),frame='fk5',unit=(u.deg,u.deg))
        
        for our_id, our_ra, our_dec, our_per in zip(our_ids, 
                                                    our_ras, 
                                                    our_decs, 
                                                    our_pers):
                                                    
            c1 = SkyCoord(str(our_ra),str(our_dec),frame='fk5',unit=(u.deg,u.deg))
            
            sep = c1.separation(c2).arcsecond
            
            if sep <= 2: 
                shared_ids = np.append(shared_ids, our_id)
                shared_our_pers = np.append(shared_our_pers, our_per)
                shared_fro_pers = np.append(shared_fro_pers, fro_per)
                id_index = np.where(our_ids==our_id)
                shared_classes = np.append(shared_classes, classes[id_index])
                break 
            
    zipped  = list(zip(shared_ids, shared_our_pers, shared_fro_pers, shared_classes))
    out_df = pd.DataFrame(zipped, columns=['ID', 'Our Per', 'Fro Per', 'qm_class'])
    out_df.to_csv(outfile, index=False)

def redo_bhard_match(our_results, bhard_results, outfile):
     # Import(s)
    import numpy as np
    import pandas as pd
    import re
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    # Action
    
    bhard_df = pd.read_csv(bhard_results)
    bhard_ras = np.array(bhard_df['RA'])
    bhard_decs = np.array(bhard_df['DEC'])
    bhard_pers = np.array(bhard_df['PER'])

    results_df = pd.read_csv(our_results)
    classes = np.array(results_df['primary_class'])
    
    our_ids = np.array(results_df['ID'])
    our_ras = np.array(results_df['RA'])
    our_decs = np.array(results_df['DEC'])
    our_pers = np.array(results_df['period'])
    
    shared_ids = np.array([])
    shared_our_pers = np.array([])
    shared_bhard_pers = np.array([])
    shared_classes = np.array([])

    counter = 0 
    
    for bhard_ra, bhard_dec, bhard_per in zip(bhard_ras, 
                                        bhard_decs, 
                                        bhard_pers):
        num = 0
        
        c2 = SkyCoord(str(bhard_ra),str(bhard_dec),frame='fk5',unit=(u.deg,u.deg))
        
        for our_id, our_ra, our_dec, our_per in zip(our_ids, 
                                                    our_ras, 
                                                    our_decs, 
                                                    our_pers):
                                                    
            c1 = SkyCoord(str(our_ra),str(our_dec),frame='fk5',unit=(u.deg,u.deg))
            
            sep = c1.separation(c2).arcsecond
            
            if sep <= 2: 
                shared_ids = np.append(shared_ids, our_id)
                shared_our_pers = np.append(shared_our_pers, our_per)
                shared_bhard_pers = np.append(shared_bhard_pers, bhard_per)
                id_index = np.where(our_ids==our_id)
                shared_classes = np.append(shared_classes, classes[id_index])
                num+=1
                counter+=1 
        
        if num > 1: 
            print(num)
            
    print(counter)
            
        
            
    zipped  = list(zip(shared_ids, shared_our_pers, shared_bhard_pers, shared_classes))
    out_df = pd.DataFrame(zipped, columns=['ID', 'Our Per', 'Bhard Per', 'qm_class'])
    out_df.to_csv(outfile, index=False)
    print(zipped)
    print(out_df)
        
def plot_comparisons(fro_csv, bhard_csv):
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    
    # Action
    fro_df = pd.read_csv(fro_csv)
    our_pers_fro = np.array(fro_df['Our Per'])
    fro_pers = np.array(fro_df['Fro Per'])
    fro_classes = np.array(fro_df['qm_class'])
    
    bhard_df = pd.read_csv(bhard_csv)
    our_ids_bhard = np.array(bhard_df['ID'])
    our_pers_bhard = np.array(bhard_df['Our Per'])
    bhard_pers = np.array(bhard_df['Bhard Per'])
    bhard_classes = np.array(bhard_df['qm_class'])
    
    for our_per, bhard_per, shared_id in zip(our_pers_bhard, bhard_pers, our_ids_bhard):
        if our_per > bhard_per*1.05 or our_per < bhard_per*0.95: 
            print('bhard missed')
            
    for our_per, fro_per in zip(our_pers_fro, fro_pers):
        if our_per > fro_per*1.05 or our_per < fro_per*0.95: 
            print('fro missed')
    
    
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
    
    
    # actual scatter points
    # Periods
    fro_mask = np.where(fro_classes=='p')
    bhard_mask = np.where(bhard_classes=='p') 
    
    plt.scatter(our_pers_fro[fro_mask], fro_pers[fro_mask], color='lightsteelblue', edgecolors='black', linewidths=0.5, label='Froebrich et al.', zorder=2.5)
    plt.scatter(our_pers_bhard[bhard_mask], bhard_pers[bhard_mask], color='#408ee0', edgecolors='black', linewidths=0.5, label='Bhardwaj et al.', zorder=2.5)
    
    # Timescales
    fro_mask = np.where(fro_classes!='p')
    bhard_mask = np.where(bhard_classes!='p') 
    
    plt.scatter(our_pers_fro[fro_mask], fro_pers[fro_mask], color='lightsteelblue', edgecolors='black', linewidths=0.5, marker='s', zorder=2.5)
    plt.scatter(our_pers_bhard[bhard_mask], bhard_pers[bhard_mask], color='#408ee0', edgecolors='black', linewidths=0.5, marker='s', zorder=2.5)
    
    # rest of plot
    
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
    #plt.show()
    plt.savefig('./referee_changes/redo_period_recovery/period_comparison.png', dpi=200, bbox_inches='tight')
    
our_results = './referee_changes/new_names/renamed.csv'
froebrich_raw = './recovered/2.0/efforts/august21/aug23-29/froebrich/froebrich_raw_data.txt'
outfile = './referee_changes/redo_period_recovery/bhard_matched.csv'
bhard_csv = './recovered/2.0/efforts/kappa/Bhardwaj_comparison/Bhardwaj_results.csv'


matched_csv = './referee_changes/redo_period_recovery/matched.csv'

match_to_file(our_results, froebrich_raw, matched_csv)
redo_bhard_match(our_results, bhard_csv, outfile)

plot_comparisons(matched_csv, bhard_csv)
