from types import new_class


def calculate_new_m(mags):
        # Import(s)
        import numpy as np

        # Action
        (tenth, ninetieth) = np.percentile(mags, [10, 90])
        mean = np.mean(mags[np.logical_or(mags>ninetieth, mags<tenth)])

        m = (mean-np.median(mags))/np.std(mags)
        return m
        
def prepare_for_vetting(results_file, data_dir, plots_dir):
    # Import(s)
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Action
    
    results_df = pd.read_csv(results_file)
    ids = np.array(results_df['ID'])
    old_ms = np.array(results_df['M'])
    qs = np.array(results_df['Q'])
    old_primary_classes = np.array(results_df['primary_class'])
    old_secondary_classes = np.array(results_df['secondary_class'])
    
    new_ms = np.array([])
    
    def make_vetting_plot(id, mjd, mag, q, old_m, new_m, old_primary_class, new_primary_class): 
        
        plt.rcParams['font.family']='serif'
        
        plt.scatter(mjd, mag, s=4, color='blue')
        median_y = np.median(mag)
        plt.axhline(y=median_y, color='green', ls='--', lw=3, label='Median')
        mean_y = np.mean(mag)
        plt.axhline(y=mean_y, color='red', ls='--', lw=3, label='Mean')
        
        title_str = id+'; Q: ' + str(round(q, 2)) + '\nOld M: '+str(round(old_m, 3))
        title_str += '; New M: '+str(round(new_m, 3)) + '\n'
        title_str += 'Old Primary: '+old_primary_class.upper() + '; Old Secondary: '
        if type(old_secondary_class)==str: 
            title_str += old_secondary_class.upper()
        
        plt.gca().set(xlabel='Time (MJD)', ylabel='Mag (r)', title=title_str)
        plt.gca().invert_yaxis()
        plt.legend()
        
        #plt.show()
        plt.savefig(os.path.join(plots_dir, id+'.png'), bbox_inches='tight')
        plt.clf()
        plt.close()
        
    for id, old_m, q, old_primary_class, old_secondary_class in zip(ids, old_ms, qs, old_primary_classes, old_secondary_classes):
        data_df = pd.read_csv(os.path.join(data_dir, id+'_r.csv'))
        mjds = np.array(data_df['mjd'])
        mags = np.array(data_df['mag'])
        
        new_m = calculate_new_m(mags)
        
        new_ms = np.append(new_ms, new_m)
        
        # Make plot
        '''
        if -0.25 < new_m < 0.25: 
            if old_m > 0.25 or old_m < -0.25: 
                make_vetting_plot(id, mjds, mags, q, old_m, new_m, old_primary_class, old_secondary_class)
        
        elif new_m > 0.25 or new_m < -0.25: 
            if -0.25 < old_m < 0.25: 
                make_vetting_plot(id, mjds, mags, q, old_m, new_m, old_primary_class, old_secondary_class)
        '''
        
        
    out_df = pd.DataFrame(list(zip(ids, new_ms)), columns=['ID','new_m'])
    out_df.to_csv('new_ms.csv', index=False)
        
        
    # within iterations 
    
    # first add to df with columns: ID, New_M
    
    # Then make plot for vetting for objects outside ranges. 
'''
results_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv'
data_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/data/AUGUST_5th/light_curves'
plots_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/september21/fixing_m/plots_for_vetting'

prepare_for_vetting(results_file, data_dir, plots_dir)
'''

def integrate_new_classes(old_results_file, new_results_file, new_m_values, new_classifications):
    # Import(s)
    import pandas as pd 
    import numpy as np
    
    # Action
    old_df = pd.read_csv(old_results_file)
    
    new_m_df = pd.read_csv(new_m_values)
    new_classes_df = pd.read_csv(new_classifications)    
    
    new_df = pd.merge(old_df, new_m_df, on='ID')
    
    new_df = new_df.drop('M', 1)
    new_df = new_df.rename(columns={"new_m":"M"})
    
    # get new qm classes mixed in
    
    ids = np.array(new_df['ID'])
    reclassified_ids = np.array(new_classes_df['ID'])
    new_classes = np.array(new_classes_df['class'])
    
    old_primary_classes = np.array(new_df['primary_class'])    
    old_secondary_classes = np.array(new_df['secondary_class'])
    new_primary_classes = np.array([])
    new_secondary_classes = np.array([])
    
    for id, old_primary, old_secondary in zip(ids, old_primary_classes, old_secondary_classes):
        if id in reclassified_ids:
            for r_id, new_class in zip(reclassified_ids, new_classes):
                if id == r_id: 
                    if type(new_class)==str:
                        if ':' in new_class: 
                            new_list = new_class.split(':')
                            new_primary_classes = np.append(new_primary_classes, new_list[0])
                            new_secondary_classes = np.append(new_secondary_classes, new_list[1])
                        
                        else: 
                            
                            if new_class == '0': 
                                new_primary_classes = np.append(new_primary_classes, old_primary)
                                new_secondary_classes = np.append(new_secondary_classes, old_secondary)
                            else: 
                                new_primary_classes = np.append(new_primary_classes, new_class)
                                new_secondary_classes = np.append(new_secondary_classes, pd.NA)
                  
        else:             
            new_primary_classes = np.append(new_primary_classes, old_primary)
            new_secondary_classes = np.append(new_secondary_classes, old_secondary)
    
    new_df = new_df.drop(['primary_class', 'secondary_class'], axis=1)
    
    new_df['primary_class'] = new_primary_classes
    new_df['secondary_class'] = new_secondary_classes
    
    new_df.to_csv(new_results_file, index=False)
    
    
new_ms = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/september21/fixing_m/new_ms.csv'
old_df = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/august21/AUGUST_FINAL_RESULTS/merged.csv'
final_df = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/september21/SEPTEMBER_FINAL_RESULTS/merged_results.csv'
new_classes_df = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/september21/fixing_m/vetting_results.csv'

integrate_new_classes(old_df, final_df, new_ms, new_classes_df)