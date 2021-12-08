def test_plot():
    #Import(s)
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    #Action

    #Plot
    plt.rcParams['font.family'] = 'serif'
    left, width = 0.25, 0.5
    bottom, height = 0.25, 0.5
    right = left+width
    top = bottom + height
    
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])

    p = patches.Rectangle((left,bottom),width,height,fill=False,transform=ax.transAxes,clip_on=False)
    ax.add_patch(p)

    ax.text(right,top,'Aperiodic',horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    ax.text(right,top,'Quasi-Periodic',horizontalalignment='center',verticalalignment='bottom',transform=ax.transAxes)
    ax.text(left,top,'Periodic',horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes)

    ax.set_xlabel('Quasi-Periodicity (Q)')
    
    plt.show()



def test_plot_2(data_file):
    #Import(s)
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    from shutil import copyfile
    from matplotlib.ticker import MultipleLocator

    #Action

    df = pd.read_csv(data_file)
    ids = list(df['ID'])
    q = list(df['Q'])
    m = list(df['M'])
    '''
    for id, quasi, asym in zip(ids,q,m):
        if 0.15 < quasi < 0.85:
            if asym < -0.25:
                copyfile('/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/plots/+++.png'.replace('+++',id),'/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/q_vs_m_results/quasi_periodic_bursters/+++.png'.replace('+++',id))
            elif asym > 0.25: 
                copyfile('/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/plots/+++.png'.replace('+++',id),'/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/q_vs_m_results/quasi_periodic_dippers/+++.png'.replace('+++',id))
    '''
            

    #Plot
    #sns.set_style('white')
    plt.rcParams['font.family'] = 'serif'

    plt.scatter(q,m,marker='.',color='slategray',s=2)

    plt.figtext(0.18,0.9,'Periodic',ha='center',fontsize='small')
    plt.figtext(0.5,0.9,'Quasi-Periodic',ha='center',fontsize='small')
    plt.figtext(0.845,0.9,'Aperiodic',ha='center',fontsize='small')

    plt.figtext(0.91,0.685,'Bursting',rotation=270,fontsize='small')
    plt.figtext(0.91,0.43,'Symmetric',rotation=270,fontsize='small')
    plt.figtext(0.91,0.22,'Dipping',rotation=270,fontsize='small')

    plt.axhline(y=0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axhline(y=-0.25,xmin=0,xmax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=0.15,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.axvline(x=0.85,ymin=-1,ymax=1,linewidth=0.8,color='black',linestyle='--')
    plt.xlim(0,1)
    plt.ylim(-1,1)
    plt.gca().invert_yaxis()
    plt.xlabel('Quasi-Periodicity (Q)')
    plt.ylabel('Flux Asymmetry (M)')
    plt.axes().yaxis.set_minor_locator(MultipleLocator(0.05))
    plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
    plt.show()

def results_table(key_file,qm_data_file,odr_data_file,ls_data_file):
    #Import(s)
    import pandas as pd
    from progress.bar import Bar

    #Action
    key_df = pd.read_csv(key_file)
    key_ids = list(key_df['ID'])
    ras = list(key_df['RA'])
    decs = list(key_df['DEC'])

    qm_df = pd.read_csv(qm_data_file)
    qm_ids = list(qm_df['ID'])
    q = list(qm_df['Q'])
    m = list(qm_df['M'])

    ls_df = pd.read_csv(ls_data_file)
    ls_ids = list(ls_df['ID'])
    periods = list(ls_df['Period'])
    powers = list(ls_df['Power'])
    faps = list(ls_df['99.9% FAP'])

    odr_df = pd.read_csv(odr_data_file)
    odr_ids = list(odr_df['star Identifier'])
    slopes = list(odr_df['odr slope'])
    slope_errors = list(odr_df['odr slope error'])

    final_ids = []
    final_ras = []
    final_decs = []
    final_types = []
    final_periods = []
    final_slopes = []

    bar = Bar('Processing...',max=len(qm_ids))
    
    for qm_id in qm_ids:
        if qm_id in ls_ids and qm_id in odr_ids:
            final_ids.append(qm_id)
            qm_id_idx = qm_ids.index(qm_id)
            q_val = q[qm_id_idx]
            m_val = m[qm_id_idx]

            if m_val < -0.25:
                final_types.append('B')
            elif q_val < 0.15 and -0.25 < m_val < 0.25:
                final_types.append('P')
            elif 0.15 < q_val < 0.85 and -0.25 < m_val < 0.25:
                final_types.append('QPS')
            elif q_val > 0.85 and -0.25 < m_val < 0.25:
                final_types.append('S')
            elif 0.15 < q_val < 0.85 and m_val > 0.25:
                final_types.append('QPD')
            elif q_val > 0.85 and m_val > 0.25:
                final_types.append('APD')
            else: 
                final_types.append('U')
            
            ls_id_idx = ls_ids.index(qm_id)
            period = periods[ls_id_idx]
            fap = faps[ls_id_idx]
            power = powers[ls_id_idx]

            if power > fap:
                final_periods.append(period)
            else:
                final_periods.append('NA')
            
            odr_id_idx = odr_ids.index(qm_id)
            odr_slope = slopes[odr_id_idx]
            odr_error = slope_errors[odr_id_idx]

            if odr_slope != 0: 
                if 10 > 100*(odr_error/odr_slope):
                    final_slopes.append(odr_slope)
                else: 
                    final_slopes.append('NA')
            else: 
                final_slopes.append('NA')

            key_id_idx = key_ids.index(qm_id)
            final_ras.append(ras[key_id_idx])
            final_decs.append(decs[key_id_idx])
            bar.next()
    
    zipped = list(zip(final_ids,final_ras,final_decs,final_types,final_periods,final_slopes))
    out_df = pd.DataFrame(data=zipped,columns=['ID','RA','DEC','TYPE','PERIOD','SLOPE'])
    out_df.to_csv('/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/results_table/results_table.csv',index=False)
    bar.finish()

def latex_the_table(raw_table,out_file):
    #Import(s)
    import pandas as pd

    #Action
    df = pd.read_csv(raw_table)
    ids = list(df['ID'])
    ras = list(df['RA'])
    decs = list(df['DEC'])
    types = list(df['TYPE'])
    periods = list(df['PERIOD'])
    slopes = list(df['SLOPE'])

    with open(out_file,'w') as f:
        for id, ra, dec, type, period, slope in zip(ids, ras, decs, types, periods, slopes): 
            if slope > 399:
                slope = 'NA'
            out_line = str(id.replace('_',' ')) + ' & ' + str(type) + ' & ' + str(round(period,3)) + ' & ' + str(round(slope,3)) 
            out_line = out_line + r' \\'
            f.write(out_line+'\n')

'''
key = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/key.csv'
qm_results = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/q_vs_m_results.csv'
odr_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/535_odr_results.csv'
ls_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/ls_results.csv'

results_table(key_file=key,qm_data_file=qm_results,odr_data_file=odr_file,ls_data_file=ls_file)
'''

latex_the_table(raw_table='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/results_table/results_table.csv',out_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/beta/results_table/latex_results.txt')