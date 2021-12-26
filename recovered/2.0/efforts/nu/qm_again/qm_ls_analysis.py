def evaluate_qm_classes(results_file, classes_file, merged_file, latex_file): 
    # Import(s)
    import pandas as pd
    import numpy as np
    
    # Definitions
    def make_plot(id, per, q, m, qm_class, 
                  data_temp, plot_dir): 
        
        # Import(s)
        import os
        import matplotlib.pyplot as plt
        from matplotlib.ticker import AutoMinorLocator, MultipleLocator
        
        # Action
        plt.rcParams['font.family']='serif'
        plt.rcParams['font.size']=5
        plt.rcParams['figure.figsize']=(3, 1.5)
        fig, ax = plt.subplots(1, 1)
        
        title_str = id.replace('_',' ')+' '+ ' Q='+ str(round(q, 2)) 
        title_str = title_str +'; M='+str(round(m, 2)) 
        title_str = title_str + ' ('+str(qm_class.upper())+')'
        fig.suptitle(title_str)
        
        ## Light curve
        
        data_df = pd.read_csv(data_temp.replace('+++', id))
        mjds = np.array(data_df['mjd'])
        mags = np.array(data_df['mag'])
        magerrs = np.array(data_df['magerr'])
        
        ### Fold if periodic 
        if qm_class == 'p':
            phased_dates = (np.mod(mjds, per))/per 
            phased_dates_cycle_2 = phased_dates + 1
            ax.scatter(phased_dates, mags, color='#408ee0', s=3)
            ax.scatter(phased_dates_cycle_2, mags, color='#408ee0', s=3)
            ax.set_xlabel('Phase (P='+str(round(per,2))+' d)')
            ax.set_ylabel('r (mag)')
            ax.xaxis.set_major_locator(MultipleLocator(0.25))
            
        else: 
            ax.errorbar(mjds, mags, yerr=magerrs, lw=0,elinewidth=0.5, color = "#408ee0")
            ax.scatter(mjds,mags,s=2, color = "#408ee0")
            ax.set_ylabel('r (mag)')
            ax.set_xlabel('Time (MJD)')
            #ax.locator_params(axis='y', nbins=5)
            ax.xaxis.set_major_locator(MultipleLocator(200))
        
        ax.xaxis.set_minor_locator(AutoMinorLocator())    
        ax.tick_params(axis='both',which='both',direction='in')
        ax.tick_params(which='both',bottom=True,top=True,left=True,right=True)
        ax.tick_params(labelbottom=True,labeltop=False,labelleft=True,labelright=False)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.invert_yaxis()
        
        plot_path = os.path.join(plot_dir, id)    
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close()
        

        
    # Action
    classes_df = pd.read_csv(classes_file) 
    results_df = pd.read_csv(results_file)
    
    merged_df = pd.merge(left=results_df, right=classes_df, sort=True, on='ID')
    
    for id, per, q, m, qm_class in zip(merged_df['ID'], merged_df['PER'],
                                       merged_df['Q'], merged_df['M'], 
                                       merged_df['primary_class']):
        
        data_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/APRIL_12th/light_curves/+++_r.csv'
        plot_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/routine_results/appendix_plots'
    
        make_plot(id, per, q, m, qm_class, data_temp, plot_dir)
        
    ### SAVE MERGED RESULTS
    merged_df.to_csv(merged_file, index=False)
    
    ### FORMAT FOR LaTeX Results Table 
    latex_df = merged_df
    latex_df.RA = latex_df.RA.round(5)
    latex_df.DEC = latex_df.DEC.round(5)
    latex_df.PER = latex_df.PER.round(3)
    latex_df.Q = latex_df.Q.round(2)
    latex_df.M = latex_df.M.round(2)
    latex_df.NU = latex_df.NU.round(3)
    latex_df.rMean = latex_df.rMean.round(2)
    latex_df.primary_class = latex_df.primary_class.str.upper()
    latex_df.secondary_class = latex_df.secondary_class.str.upper()
    latex_df = latex_df.replace(np.nan, '-')
    
    file_lines = []
    
    # Write result strings to text file
    for counter in range(len(latex_df.ID)):
        result_string = ''
        df_columns = list(latex_df.columns)
        for col_name in df_columns:
            result_string = result_string+str((latex_df[col_name])[counter])
            if df_columns.index(col_name) < len(df_columns)-1:
                result_string = result_string + ' & '
        
        result_string = result_string+r'\\' # Neeed the double bar
        result_string = result_string.replace('_',' ')
        
        file_lines.append(result_string)
    
    with open(latex_file, 'w') as f:
        f.write(' \n'.join(file_lines))
        
            
                    
results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/routine_results/results.csv'
classes_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/routine_results/classes.csv'
merged_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/routine_results/finalized_results/merged_results.csv'
latex_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/routine_results/finalized_results/table_formatted.txt'
evaluate_qm_classes(results_file, classes_file, merged_file, latex_file)


# plots: /home/thaddaeus/FMU/HRL/LAH2.0/efforts/nu/qm_again/routine_results/classes.csv