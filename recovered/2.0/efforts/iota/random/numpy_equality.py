# Import(s)
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Action
size = 1000 
qs = np.random.randint(0,100,size=size)
even_qs_mask = np.where(qs%2==0)

ms = np.random.randint(0,25,size=size)
odd_ms_mask = np.where(ms%2!=0)

intersection_mask = np.intersect1d(odd_ms_mask,even_qs_mask)
sns.jointplot(qs[intersection_mask],ms[intersection_mask])
plt.show()