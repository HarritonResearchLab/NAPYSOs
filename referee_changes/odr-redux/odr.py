import numpy as np
import pandas as pd
import re
import csv
import glob
import os
from pathlib import Path
from scipy import odr


def find_odr(cmd_data):
    # get cmd data
    df = pd.read_csv(cmd_data)

    # get g and g-r/color data
    g_minus_r = np.array(df['color'])
    g_minus_r_magerrs = np.array(df['color_err'])
    g_vals = np.array(df['g_mag'])
    cmd_g_magerrs = np.array(df['g_magerr'])

    # set up for odr
    # set arrays to lists
    g_vals = list(g_vals)
    g_minus_r_magerrs = list(g_minus_r_magerrs)
    cmd_g_magerrs = list(cmd_g_magerrs)
    g_minus_r = list(g_minus_r)

    # exclude < 2.5th and > 97.5th percentile values
    # find percentiles
    g_percentiles = np.percentile(g_vals, q=[2.5, 97.5])
    g_minus_r_percentiles = np.percentile(g_minus_r, [2.5, 97.5])

    # remove the outside 5 percent
    for item in g_vals:
        if item < g_percentiles[0] or item > g_percentiles[1]:
            item_index = g_vals.index(item)
            g_vals.remove(item)
            g_minus_r.remove(g_minus_r[item_index])
            cmd_g_magerrs.remove(cmd_g_magerrs[item_index])
            g_minus_r_magerrs.remove(g_minus_r_magerrs[item_index])

    for item in g_minus_r:
        if item < g_minus_r_percentiles[0] or item > g_minus_r_percentiles[1]:
            item_index = g_minus_r.index(item)
            g_minus_r.remove(item)
            g_vals.remove(g_vals[item_index])
            cmd_g_magerrs.remove(cmd_g_magerrs[item_index])
            g_minus_r_magerrs.remove(g_minus_r_magerrs[item_index])

    # return to arrays
    g_vals = np.array(g_vals)
    g_minus_r_magerrs = np.array(g_minus_r_magerrs)
    cmd_g_magerrs = np.array(cmd_g_magerrs)
    g_minus_r = np.array(g_minus_r)

    # odr
    def odr_line(B, x):
        y = B[0] * x + B[1]
        return y

    def perform_odr(x, y, xerr, yerr):
        linear = odr.Model(odr_line)
        mydata = odr.Data(x, y, wd=1. / xerr, we=1. / yerr)
        myodr = odr.ODR(mydata, linear, beta0=[0, 0])
        output = myodr.run()
        return output

    m = np.array([])
    b = np.array([])
    for i in np.arange(0, 1000):
        indices = np.arange(0, len(g_vals))
        np.random.shuffle(indices)
        ind = indices[0:int(round(len(g_vals) / 2, 0))]  # dividing by integer on purpose.
        out = perform_odr(g_minus_r[ind], g_vals[ind], g_minus_r_magerrs[ind], cmd_g_magerrs[ind])
        m = np.append(m, out.beta[0])
        b = np.append(b, out.beta[1])


    a = np.arctan(m)  # in radian
    # get median and symmetric 68% interval for m, b and alpha:
    m_median = np.median(m)
    m_down = np.sort(m)[int(round(0.16 * len(m)))]
    m_up = np.sort(m)[int(round(0.84 * len(m)))]
    m_error = np.mean([abs(m_down - m_median), abs(m_up - m_median)])
    # print (m_median, m_up, m_down, m_error)
    b_median = np.median(b)
    b_down = np.sort(b)[int(round(0.16 * len(b)))]
    b_up = np.sort(b)[int(round(0.84 * len(b)))]
    b_error = np.mean([abs(b_down - b_median), abs(b_up - b_median)])
    # print (b_median, b_up, b_down, b_error)
    a_median = np.median(a)
    a_down = np.sort(a)[int(round(0.16 * len(a)))]
    a_up = np.sort(a)[int(round(0.84 * len(a)))]
    a_error = np.mean([abs(a_down - a_median), abs(a_up - a_median)])

    pythag = np.sqrt((max(g_vals)-min(g_vals))*(max(g_vals)-min(g_vals))+(max(g_minus_r)-min(g_minus_r))*(max(g_minus_r)-min(g_minus_r)))

    plot_title = re.search(r'[^\/\\]+(?=_cmd.csv)', cmd_data)   # this regex matches for starname_cmd.csv, use [^\/\\]+(?=.csv) for starname.csv

    return plot_title.group(), m_median, m_error, b_median, b_error, np.rad2deg(a_median), np.rad2deg(a_error), pythag


base_path = r'C:\Users\sunto\PycharmProjects\pythonProject1\referee_cmd_files'

with open(os.path.join(base_path, "referee_odr.csv"), "r") as f:
    d = f.read()
with open(os.path.join(base_path, "referee_odr.csv"), "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["star Identifier", "odr slope", "odr slope error", "intercept", "intercept error", "angle", "angle error", "pythag"])
    for cmd_data in glob.glob(os.path.join(base_path, "*_cmd.csv")):
        cmd_data_f_n = Path(cmd_data).stem
        print(cmd_data_f_n)
        name, slope, slope_error, intercept, intercept_error, angle, angle_error, pythag = find_odr(cmd_data)
        writer.writerow([cmd_data_f_n, slope, slope_error, intercept, intercept_error, angle, angle_error, pythag])
