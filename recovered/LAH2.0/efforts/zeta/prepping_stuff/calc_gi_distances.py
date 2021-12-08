def calc_gi_dists(key, r_path, g_path, placeholder):
    # Import(s)
    import pandas as pd
    import numpy as np
    from astropy.io import ascii
    from personalastropy.ysospy.interpolation import returnGoodRegions
    from personalastropy.ysospy.handy_scripts import sortData
    import matplotlib.pyplot as plt
    from progress.bar import Bar
    import seaborn as sns

    # Action

    df = pd.read_csv(key)

    seps = []  # list of all the observation separations

    ids = list(df['ID'])
    bar = Bar('Processing', max=len(ids))
    for id in ids:
        r_file = r_path.replace(placeholder, id)
        r = ascii.read(r_file, format='ipac', delimiter='|')
        r_mjds = list(r['mjd'])
        r_mags = list(r['mag'])
        r_magerrs = list(r['magerr'])

        g_file = g_path.replace(placeholder, id)
        g = ascii.read(g_file, format='ipac', delimiter='|')
        g_mjds = list(g['mjd'])
        g_mags = list(g['mag'])
        g_magerrs = list(g['magerr'])

        cleaned_r_mjds = []
        cleaned_r_mags = []
        cleaned_r_magerrs = []

        for item in r_mjds:
            if item > 58456 or item < 58448:
                item_index = r_mjds.index(item)
                cleaned_r_mjds.append(item)
                cleaned_r_mags.append(r_mags[item_index])
                cleaned_r_magerrs.append(r_magerrs[item_index])

        red_sorter = sortData(x=cleaned_r_mjds, y=[
                              cleaned_r_mags, cleaned_r_magerrs])

        srd = red_sorter[0]
        srm = red_sorter[1][0]
        sre = red_sorter[1][1]

        green_sorter = sortData(x=g_mjds, y=[g_mags, g_magerrs])

        sgd = green_sorter[0]
        sgm = green_sorter[1][0]
        sge = green_sorter[1][1]

        # Process and save cmd data to csv
        # good green intervals func
        gis = returnGoodRegions(x=sgd, y=sgm, max_sep=3, min_card=2)

        green_date_intervals = gis[1]
        for i in green_date_intervals:
            #i = list(i)
            max_index = len(i)-2
            for date_indice, date in enumerate(list(i)):
                if date_indice < max_index:
                    next_date_index = date_indice+1
                    sep = (i[next_date_index]-date)
                    seps.append(sep)

        bar.next()

    with open('seps.txt', 'w') as f:
        for i in seps:
            f.write(str(i)+'\n')

    # Some calculations
    same_or_below_indices = np.where(np.array(seps) <= 2)
    percentile = len(same_or_below_indices)/len(seps)

    # Plot
    sns.set_style('darkgrid')
    plt.rcParams['font.family'] = 'serif'
    plt.hist(x=seps, bins=150, label='Max sep: 3d')

    plt.axvline(x=2, linestyle='--', label=str(percentile))
    # sns.displot(seps,kind='kde')
    plt.title('Justification for maximum interpolation distance')
    plt.xlabel('Observation separation')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    plt.show()


'''
r_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++++_r.tbl'
g_temp = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/data/lightcurves/+++++_g.tbl'
key_path = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/for_distribution/key.csv'

calc_gi_dists(key=key_path,r_path=r_temp,g_path=g_temp,placeholder='+++++')
'''


def just_analyze(seps_file):
    # Import(s)
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    # Action
    seps = []
    with open(seps_file, 'r') as f:
        for line in f:
            line = line.replace('\n', '')
            seps.append(float(line))

    # Some calculations
    same_or_below1 = []
    same_or_below2 = []
    for i in seps:
        if i <= 2.25:
            same_or_below1.append(i)
            if i <= 2.0:
                same_or_below2.append(i)
    percentile_1 = len(same_or_below1)/len(seps)*100
    percentile_2 = len(same_or_below2)/len(seps)*100

    # Plot
    sns.set_style('darkgrid')
    plt.rcParams['font.family'] = 'serif'
    plt.hist(x=seps, bins=150, label='Max Sep (d): 3')

    plt.axvline(x=2.25, linestyle='--', lw=0.5,
                label=str(round(percentile_1, 1))+'th percentile')
    plt.axvline(x=2, linestyle='--', lw=0.5,
                label=str(round(percentile_2, 1))+'th percentile')
    plt.title('Justification for maximum interpolation distance')
    plt.xlabel('Observation separation')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right', fontsize='small')
    plt.show()


just_analyze(seps_file='seps.txt')
