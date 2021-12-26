from sys import base_prefix
import numpy as np
from numpy.core.numeric import indices

ls_freqs = np.linspace(0,15,15)
bap = np.array([1,3,5,6,9,12])
dv = 0.9
indices_to_remove = []
for f in bap: 
    mask = np.invert(np.logical_and((ls_freqs+dv)>f, (ls_freqs-dv)<f))
    idx = np.argwhere(mask==False)
    for i in idx:
        indices_to_remove.append(int(i))


indices_to_remove = np.array(sorted(indices_to_remove,reverse=True))

cleaned_ls_freqs = np.delete(ls_freqs,indices_to_remove)
print(cleaned_ls_freqs)

