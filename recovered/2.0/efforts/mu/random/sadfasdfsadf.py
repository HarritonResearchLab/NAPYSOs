import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

x1 = np.random.normal(0,1,20)
y1 = np.random.normal(0,1,20)

x2 = np.random.normal(3, 1.2, 20)
y2 = np.random.normal(3,1.2,20)

x3 = np.random.normal(1, 1.1, 20)
y3 = np.random.normal(1, 1.1, 20)

colors = sns.color_palette(n_colors=3, as_cmap=True)

plt.scatter(x1,y1,color=colors[0])
plt.scatter(x2,y2,color=colors[1])
plt.scatter(x3,y3,color=colors[2])
plt.show()