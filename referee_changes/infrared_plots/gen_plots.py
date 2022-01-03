"""Modified from snippet used for Bredall et al. (2020) Fig. 7"""
# Import(s)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

og_x_g = ['Kmag', 'W3mag']  # set to ZEROS if you want zeros
og_y_g = ['Kmag', 'W4mag']  # set to ZEROS if you want zeros

color_marker_dict = {
    'p': ['C0', 's'],
    'qps': ['C1', 'o'],
    's': ['C2', '*'],
    'apd': ['C3', 'D'],
    'qpd': ['C4', 'h'],
    'b': ['C5', 'd'],
    'l': ['black', 'p'],
    'mp': ['C6', 'P']
}


def gen_infared_plot(label_x: str,
                     label_y: str,
                     of: str,
                     og_x: list = og_x_g,
                     og_y: list = og_y_g):
    i_df = pd.read_csv('updated_object_infrared.csv')
    i_df.set_index('ID_1')
    i_df['ZEROS'] = np.float64(0.0)

    for x in set(['ID_1', 'primary_class'] + og_x + og_y):
        i_df = i_df[i_df[x].notna()]
        if x in ['ID_1', 'primary_class', 'ZEROS']:
            continue
        i_df[x] = i_df[x].astype(np.float64)

    i_df = i_df[i_df['primary_class'] != 'u']

    p = i_df[i_df['primary_class'] == 'p']
    q = i_df[i_df['primary_class'] == 'q']
    qps = i_df[i_df['primary_class'] == 'qps']
    qpd = i_df[i_df['primary_class'] == 'qpd']
    b = i_df[i_df['primary_class'] == 'b']
    s = i_df[i_df['primary_class'] == 's']
    l = i_df[i_df['primary_class'] == 'l']
    apd = i_df[i_df['primary_class'] == 'apd']
    mp = i_df[i_df['primary_class'] == 'mp']

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    lw = 0.5
    plt.scatter(p[og_x[0]] - p[og_x[1]],
                p[og_y[0]] - p[og_y[1]],
                label='Periodic',
                linewidths=lw,
                color=color_marker_dict['p'][0],
                marker=color_marker_dict['p'][1],
                edgecolors='black')
    plt.scatter(qps[og_x[0]] - qps[og_x[1]],
                qps[og_y[0]] - qps[og_y[1]],
                label='Quasi-periodic',
                linewidths=lw,
                color=color_marker_dict['qps'][0],
                marker=color_marker_dict['qps'][1],
                edgecolors='black')
    plt.scatter(s[og_x[0]] - s[og_x[1]],
                s[og_y[0]] - s[og_y[1]],
                label='Stochastic',
                linewidths=lw,
                color=color_marker_dict['s'][0],
                marker=color_marker_dict['s'][1],
                edgecolors='black')
    plt.scatter(apd[og_x[0]] - apd[og_x[1]],
                apd[og_y[0]] - apd[og_y[1]],
                label='Aperiodic dipper',
                linewidths=lw,
                color=color_marker_dict['apd'][0],
                marker=color_marker_dict['apd'][1],
                edgecolors='black')
    plt.scatter(qpd[og_x[0]] - qpd[og_x[1]],
                qpd[og_y[0]] - qpd[og_y[1]],
                label='Quasi-periodic dipper',
                linewidths=lw,
                color=color_marker_dict['qpd'][0],
                marker=color_marker_dict['qpd'][1],
                edgecolors='black')
    plt.scatter(b[og_x[0]] - b[og_x[1]],
                b[og_y[0]] - b[og_y[1]],
                label='Burster',
                linewidths=lw,
                color=color_marker_dict['b'][0],
                marker=color_marker_dict['b'][1],
                edgecolors='black')
    plt.scatter(l[og_x[0]] - l[og_x[1]],
                l[og_y[0]] - l[og_y[1]],
                label='Long-timescale',
                color=color_marker_dict['l'][0],
                marker=color_marker_dict['l'][1],
                linewidths=lw,
                edgecolors='black')
    plt.scatter(mp[og_x[0]] - mp[og_x[1]],
                mp[og_y[0]] - mp[og_y[1]],
                label='Multiperiodic',
                color=color_marker_dict['mp'][0],
                marker=color_marker_dict['mp'][1],
                linewidths=lw,
                edgecolors='black')

    # Disk category classses
    '''plt.plot((6.5, 3.5), (2.7, 1.5), c="black", linestyle="--")
    plt.plot((5, 2.5), (0, 2.5), c="black", linestyle="--")
    plt.plot((2.4, 1.5), (-0.25, 1.25), c="black", linestyle="--")
    plt.plot((0.375, 0), (-0.25, 0.5), c="black", linestyle="--")
    plt.annotate("Diskless", (-0.45, -0.45), fontsize=10)
    plt.annotate("Debris", (1, -0.25), fontsize=10)
    plt.annotate("Evolved", (3, 0), fontsize=10)
    plt.annotate("Transition", (5, 1.75), fontsize=10)
    plt.annotate("Full", (3.25, 2.75), fontsize=10)'''
    plt.legend(loc='upper left',
               fontsize=7,
               fancybox=False,
               edgecolor='black',
               shadow=False,
               ncol=2)

    plt.minorticks_on()
    plt.tick_params(axis='both', which='both', direction='in')
    plt.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    plt.tick_params(labelbottom=True,
                    labeltop=False,
                    labelleft=True,
                    labelright=False)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    # plt.xlabel(r'$K_s-[4.5\mu$'+'m] (mag)')
    # plt.ylabel(r'$K_s-[3.6\mu$'+'m] (mag)')
    # plt.xlabel(r'$K_s-[12\mu$'+'m] (mag)')
    # plt.ylabel(r'$K_s-[22\mu$'+'m] (mag)')
    # plt.xlabel(r'H-K')
    # plt.ylabel(r'J-H')
    # plt.xlabel(r'$[5.8\mu \mathrm{m}]-[8.0\mu \mathrm{m}]$ (mag)')
    # plt.ylabel(r'$[3.6\mu \mathrm{m}]-[4.5\mu \mathrm{m}]$ (mag)')
    # plt.xlabel(r'$\nu$')
    # plt.ylabel(r'$[3.4\mu \mathrm{m}]-[4.6\mu \mathrm{m}]$ (mag)')

    plt.savefig(of, dpi=600)
    plt.close()


if __name__ == '__main__':
    gen_infared_plot(
        og_x=['NU', 'ZEROS'],
        og_y=['W1mag', 'W2mag'],
        label_x=r'$\nu$',
        label_y=r'$[3.4\mu \mathrm{m}]-[4.6\mu \mathrm{m}]$ (mag)',
        of='nu_vs_34-46.png')
    gen_infared_plot(
        og_x=['__5_8_', '__8_0_'],
        og_y=['__3_6_', '__4_5_'],
        label_x=r'$[5.8\mu \mathrm{m}]-[8.0\mu \mathrm{m}]$ (mag)',
        label_y=r'$[3.6\mu \mathrm{m}]-[4.5\mu \mathrm{m}]$ (mag)',
        of='58-80_vs_36-45.png')
    gen_infared_plot(og_x=['Kmag', '__4_5_'],
                     og_y=['Kmag', '__3_6_'],
                     label_x=r'$K_s-[4.5\mu\mathrm{m}]$ (mag)',
                     label_y=r'$K_s-[3.6\mu\mathrm{m}]$ (mag)',
                     of='ks-45_vs_ks-36.png')
    gen_infared_plot(og_x=['Hmag', 'Kmag'],
                     og_y=['Jmag', 'Hmag'],
                     label_x=r'$H-K_s$ (mag)',
                     label_y=r'$J-H$ (mag)',
                     of='h-k_vs_j-h.png')
    gen_infared_plot(og_x=['Kmag', 'W3mag'],
                     og_y=['Kmag', 'W4mag'],
                     label_x=r'$K_s-[12\mu\mathrm{m}]$ (mag)',
                     label_y=r'$K_s-[22\mu\mathrm{m}]$ (mag)',
                     of='ks-12_vs_ks-22.png')
