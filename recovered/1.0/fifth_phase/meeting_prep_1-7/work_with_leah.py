#Import(s)
import matplotlib.pyplot as plt
from astropy.io import ascii

coords = []
with open('/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/coordinates.txt','r') as f:
    for line in f:
        line = line.replace('\n','')
        coords.append(line)

for coord in coords: 
    path_temp = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_+++++_=.tbl'
    g_file = path_temp.replace('+++++',coord)
    g_file = g_file.replace('=','g')
    r_file = path_temp.replace('+++++',coord)
    r_file = r_file.replace('=','r')

    r=ascii.read(r_file,format='ipac',delimiter='|')
    g=ascii.read(g_file,format='ipac',delimiter='|')

    red_mags=r['mag']
    red_dates=r['hjd']
    red_magerrs=r['magerr']

    green_mags = g['mag']
    green_dates=g['hjd']
    green_magerrs=g['magerr']

    fig, axs = plt.subplots(2,1)

    axs[0].scatter(red_dates,red_mags,color='tomato')
    axs[0].invert_yaxis()

    axs[1].scatter(green_dates,green_mags,color='green')
    axs[1].invert_yaxis()

    save_path = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/plots/+++++.png'
    file_path = save_path.replace('+++++',coord)
    plt.savefig(file_path)
    plt.clf()
