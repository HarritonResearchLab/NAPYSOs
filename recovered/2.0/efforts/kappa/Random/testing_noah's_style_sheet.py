# Import(s)
import numpy as np
import matplotlib.pyplot as plt

# Action
plt.style.use('hrl')
fig, axs = plt.subplots(2,2)


x = np.linspace(0,10,3)
axs[0,0].scatter(x,x)
axs[0,0].scatter(x,2*x)
axs[0,0].scatter(x,0.5*x)
axs[0,0].set_xlabel('x')
axs[0,0].set_ylabel('y')

axs[0,1].plot(x,x)
axs[0,1].plot(x,2*x)
axs[0,1].plot(x,0.5*x)
axs[0,1].set_xlabel('x')
axs[0,1].set_ylabel('y')

normal_x = np.random.normal(0,1,100)
axs[1,0].hist(normal_x,bins=20)
axs[1,0].set_xlabel('x')
axs[1,0].set_ylabel('Frequency')

axs[1,1].axis('off')

plt.show()
