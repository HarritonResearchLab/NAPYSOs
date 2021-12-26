from astropy.units.cgs import C


def match_objects(our_results, bhardwaj_results, out_path): 
    # Import(s)
    import pandas as pd
    import numpy as np
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    # Action
    our_df = pd.read_csv(our_results)
    our_df = our_df.dropna(subset=['PER'])
    
    our_ids = np.array(our_df['ID'])
    our_ras, our_decs = np.array(our_df['RA']), np.array(our_df['DEC'])
    
    bhard_df = pd.read_csv(bhardwaj_results)
    bhard_ras, bhard_decs = np.array(bhard_df['RA']), np.array(bhard_df['DEC'])
    bhard_pers = np.array(bhard_df['PER'])
    
    matched_ids = np.array([])
    matched_bhard_pers = np.array([])

    for bhard_ra, bhard_dec, bhard_per in zip(bhard_ras, bhard_decs, bhard_pers):
        c1 = SkyCoord(str(bhard_ra),str(bhard_dec),unit=(u.deg,u.deg))
        for id, our_ra, our_dec in zip(our_ids, our_ras, our_decs):
            c2 = SkyCoord(str(our_ra),str(our_dec),unit=(u.deg,u.deg))
            sep = c1.separation(c2).arcsecond
            if sep < 0.5: 
                matched_ids = np.append(matched_ids, id)
                matched_bhard_pers = np.append(matched_bhard_pers, bhard_per)
                break

    matched_indices = np.array([])
    
    for matched_id in matched_ids: 
        for id in our_ids: 
            if matched_id == id: 
                matched_indices = np.append(matched_indices, np.where(our_ids==id)[0])
                break
    
    matched_indices = matched_indices.astype(int)
    
    our_pers = np.array(our_df['PER'])
    our_pers = our_pers[matched_indices]
    our_classes = np.array(our_df['primary_class'])
    our_classes = our_classes[matched_indices]
    
    zipped = list(zip(matched_ids, our_pers, matched_bhard_pers, our_classes))
    df_columns = ['ID','OUR_PER', 'BHARD_PER', 'primary_class']
    combined_df = pd.DataFrame(zipped, columns=df_columns)
    combined_df.to_csv(out_path, index=False)

    
our_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Misc./finalized_results/merged_results.csv'
bhardwaj_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/kappa/Bhardwaj_comparison/Bhardwaj_results.csv'
out_path = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Period Vetting/Bhardwaj/overlap.csv'
#match_objects(our_results, bhardwaj_results, out_path)


def investigate_overlap(overlap_file):
    # Import(s)
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator

    
    # Action
    df = pd.read_csv(overlap_file)
    our_ids = np.array(df['ID'])
    our_pers = np.array(df['OUR_PER'])
    bhard_pers = np.array(df['BHARD_PER'])
    qm_classes = np.array(df['primary_class'])
    
    def cody_Q(mjd, mag, magerr, timescale, sig_factor):
        # Import(s)
        from astropy.convolution import Box1DKernel, convolve
        
        # Convert to arrays
        mjd = np.array(mjd)
        mag = np.array(mag)
        magerr = np.array(magerr)
            
        # Calculate sig
        sig = sig_factor*np.mean(magerr)

        # Create the residual curve
        phase = mjd % timescale
        mag = mag[np.argsort(phase)]

        # We use three periods and extract the middle to prevent edge effects
        three_periods = np.concatenate((mag, mag, mag))
        boxcar = Box1DKernel(len(mag) // 4)
        smooth_mag = convolve(three_periods, boxcar)
        smooth_mag = smooth_mag[np.size(mag):2*np.size(mag)]

        resid_mag = mag - smooth_mag

        q = float((np.nanstd(resid_mag) ** 2 - sig ** 2)/(np.nanstd(mag) ** 2 - sig ** 2))

        return q, resid_mag
    
    def make_comparison_plot():
        plt.rcParams['font.family']='serif'
        plt.rcParams['figure.figsize']=(5, 5)
        ax = plt.gca()
        
        plt.axline(xy1=(0, 0), slope=1, color='grey', ls='--')
        plt.axline(xy1=(0, 0), slope=0.5, color='grey', ls='--')
        plt.axline(xy1=(0, 0), slope=2, color='grey', ls='--')
        
        plt.scatter(our_pers, bhard_pers, color='#408ee0', edgecolors='black', linewidths=0.5)
        plt.xlabel('Period (d; this work)')
        plt.ylabel('Period (d; Bhardwaj+)')
        
        plt.xlim(left=0.5, right=14)
        plt.ylim(bottom=0.5, top=14)
        plt.tick_params(axis='both',which='both',direction='in')
        plt.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        plt.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        ax.xaxis.set_minor_locator(AutoMinorLocator())    
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        #plt.show()
        plt.savefig('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Period Vetting/Bhardwaj/comparison.png', dpi=400, bbox_inches='tight')
    
    def make_diagnostic_plot(id, mjd, mag, our_per, bhardwaj_per, our_q, bhardwaj_q, plot_dir):
        # Import(s)
        import os 
        
        # Action
        plt.rcParams['font.family']='serif'
        plt.rcParams['figure.figsize']=(7, 4)
        fig, axs = plt.subplots(1, 2)
        
        for i in range(2):
            axs[i].tick_params(axis='both',which='both',direction='in')
            axs[i].tick_params(which='both',bottom=True,top=True,left=True,right=True)
            axs[i].tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        
        # Our period 
        phased_dates = (np.mod(mjd, our_per))/our_per 
        phased_dates_cycle_2 = phased_dates + 1
        axs[0].scatter(phased_dates, mag, color='#408ee0', s=3)
        axs[0].scatter(phased_dates_cycle_2, mag, color='#408ee0', s=3)
        axs[0].set_xlabel('Our Folded LC (P='+str(round(our_per, 2))+')')
        axs[0].set_ylabel('r (mag)')
        axs[0].set_title('Q: '+str(round(our_q, 2)), fontsize='medium')
        
        # Rebull period 
        phased_dates = (np.mod(mjd, bhardwaj_per))/bhardwaj_per 
        phased_dates_cycle_2 = phased_dates + 1
        axs[1].scatter(phased_dates, mag, color='#408ee0', s=3)
        axs[1].scatter(phased_dates_cycle_2, mag, color='#408ee0', s=3)
        axs[1].set_xlabel('Bhardwaj+ Folded LC (P='+str(round(bhardwaj_per, 2))+')')
        axs[1].set_ylabel('r (mag)')
        axs[1].set_title('Q: '+str(round(bhardwaj_q, 2)), fontsize='medium')
        
        plt.subplots_adjust(wspace=0.4, hspace=0.3)
        
        plt.savefig(os.path.join(plot_dir, (id+'.png')))
        plt.clf()
        plt.close()
    
    for id, our_per, bhard_per in zip(our_ids, our_pers, bhard_pers):
        if bhard_per > 0.5 and np.abs(our_per-bhard_per) > 0.1*(our_per):   
            data_df = pd.read_csv('/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'.replace('+++', id))
            mjds = np.array(data_df['mjd'])
            mags = np.array(data_df['mag'])
            magerrs = np.array(data_df['magerr'])
            
            our_q = cody_Q(mjds, mags, magerrs, our_per, 1.25)[0]
            bhard_q = cody_Q(mjds, mags, magerrs, bhard_per, 1.25)[0]
            
            plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Period Vetting/Bhardwaj/investigating_differences/plots'
            make_diagnostic_plot(id, mjds, mags, our_per, bhard_per, our_q, bhard_q, plot_dir)
        
        
            
    
    #make_comparison_plot()
    
    def fold_to_bhardwajs(our_ids, bhard_pers):
        #indices = np.where(qm_classes=='p')
        #bhard_pers = bhard_pers[indices]
        #our_ids = our_ids[indices]
        
        for id, per in zip(our_ids, bhard_pers): 
            data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_7th/light_curves/+++_r.csv'
            data_df = pd.read_csv(data_temp.replace("+++", id))
            
            mjds = np.array(data_df['mjd'])
            mags = np.array(data_df['mag'])
            phased_dates = (np.mod(mjds, per))/per 
            phased_dates_cycle_2 = phased_dates + 1
            plt.scatter(phased_dates, mags, c='#408ee0')
            plt.scatter(phased_dates_cycle_2, mags, c='#408ee0')
            plt.xlabel('Phase (P='+str(per)+')')
            plt.gca().invert_yaxis()
            plt.ylabel('r (mag)')
            plot_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Period Vetting/Bhardwaj/investigating_differences/forced_to_bhards/+++.png'
            plt.savefig(plot_temp.replace('+++', id))
            plt.clf()
            plt.close()
               
    #fold_to_bhardwajs(['FHK_26', 'GDR2_2163139804529410944', '2MASS_J20505543+4417461', 'FHK_82', '2MASS_J20510029+4424364'], [0.518, 10.878, 9.383, 1.162, 1.39])
                
    
overlap_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/xi/Period Vetting/Bhardwaj/overlap.csv'
investigate_overlap(overlap_file)