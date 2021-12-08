# Import(s)
import numpy as np
import matplotlib.pyplot as plt

# Action
solar_freq = 1.0
solar_per = 1.0

harmonics = np.array(5*[1])/range(1,6)
lower_harmonics = 0.98*harmonics
upper_harmonics = 1.02*harmonics
lph = 1/upper_harmonics # lph = Lower Period Harmonics (for lower bounds of harmonic period ranges)
uph = 1/lower_harmonics # uph = Upper Period Harmonics (for upper bounds of harmonic period ranges)

for i in range(5):
    print(lph[i],uph[i])

for i in range(0,5):
    print(i)
