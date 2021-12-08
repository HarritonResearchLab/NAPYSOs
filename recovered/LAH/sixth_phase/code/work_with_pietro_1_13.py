#Import(s)
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii 
import pandas as pd

#Action

coord_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/coordinates.txt'
data_temp = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/cmd_data/+++.csv'

with open(coord_file,'r') as f:
    for line in f:
        line = line.replace('\n','')
        csv_path = data_temp.replace('+++',line)
        
        

#'/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/for_email/cmd_data/+++.csv'

