'''
import re 
import os, shutil
import requests
import gzip
import astropy 
import astroquery
import sys
import statistics as stat
import os.path
import numpy as np
import pandas as pd
from astroquery.mearth import Mearth
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
'''



def querySimbad():
    #Import(s)
    import numpy as np
    from astropy import coordinates as coord
    from astropy import units as u
    from astroquery.simbad import Simbad
    from astropy.coordinates.sky_coordinate import SkyCoord
    from personalastropy.ysospy.handy_scripts import queryCoordSimbad


    coordinates = []
    with open('/home/thaddaeus/FMU/HRL/LAH/third_phase/data/obsids.txt','r') as f:
        for line in f:
            line = line.replace('\n','')
            coordinates.append(line)

    for raw_coordinate in coordinates:
        print(queryCoordSimbad(raw_coordinate,5))
        
        
    
querySimbad()

