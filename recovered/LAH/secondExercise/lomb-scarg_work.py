
from astropy.io import ascii
import numpy as np

from personalastropy.ysospy.variability import calculateFluxAsymmetry


red_file = '/home/thaddaeus/FMU/HRL/LAH/third_phase/data/lightcurves/lc_'+'20:52:53.43+44:19:36.3'+'_r.tbl'
r = ascii.read(red_file,format='ipac',delimiter='|') 
red_mags = np.array(r['mag'])

print(calculateFluxAsymmetry(x=red_mags))
