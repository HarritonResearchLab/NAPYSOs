import numpy as np
import matplotlib.pyplot as plt
from personalastropy.ysospy.plotting_funcs import plotLightCurve
from astropy.io import ascii
import math

with open('/home/thaddaeus/FMU/HRL/LAH/fifth_phase/for_email/coordinates.txt','r') as f:
    for line in f:
        line = line.replace('\n','')
        r_file = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_+++_r.tbl'.replace('+++',line)
        g_file = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_+++_g.tbl'.replace('+++',line)

        r = ascii.read(r_file,format='ipac',delimiter='|')
        g = ascii.read(g_file,format='ipac',delimiter='|')

        r_dates = np.array(r['mjd'])
        r_mags = np.array(r['mag'])
        r_magerrs = np.array(r['magerr'])

        g_dates = np.array(g['mjd'])
        g_mags = np.array(g['mag'])
        g_magerrs = np.array(g['magerr'])

        min_rd = math.floor(np.max(r_dates))
        max_rd = math.ceil(np.max(r_dates))
        r_bins = np.linspace(min_rd,max_rd,(max_rd-min_rd-1))
        print(r_bins)

        min_gd = math.floor(np.max(g_dates))
        max_gd = math.floor(np.max(r_dates))
        g_bins = max_gd - min_gd