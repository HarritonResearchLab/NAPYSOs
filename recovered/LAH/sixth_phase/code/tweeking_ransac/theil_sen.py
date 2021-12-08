#Import(s)
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import TheilSenRegressor, HuberRegressor
import pandas as pd
import seaborn as sns

#Action

#Get data
coord = '20:51:19.43+44:19:30.5'
csv_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/meeting_prep_1-7/cmd_data/+++++.csv'.replace('+++++',coord)

df = pd.read_csv(csv_file)
g_minus_rs = np.array(df['g-r'])
g = np.array(df['g'])

g_minus_rs_2d = np.atleast_2d(g_minus_rs).T
g_2d = np.atleast_2d(g).T

theil = TheilSenRegressor()
theil.fit(g_minus_rs_2d, g_2d)

predicted_x = np.arange(min(g_minus_rs),max(g_minus_rs))[:, np.newaxis]
theil_y = theil.predict(predicted_x)

#Huber
epsilon_value = 1.03
huber = HuberRegressor(epsilon=epsilon_value)
huber.fit(g_minus_rs_2d,g_2d)
predicted_x = np.arange(min(g_minus_rs),max(g_minus_rs))[:, np.newaxis]

huber_y = huber.predict(predicted_x)

#Plot

fig, axs = plt.subplots(2,1,sharex=True,sharey=True)

sns.set_style('white')

axs[0].scatter(g_minus_rs, g,marker='.',color='indianred')
axs[0].plot(predicted_x, theil_y, c='cornflowerblue',label='Theil-Sen')
axs[0].invert_yaxis()
box = axs[0].get_position()
axs[0].set_position([box.x0,box.y0,box.width*0.8,box.height])
axs[0].legend(bbox_to_anchor=(1,0.5),loc='center left',fontsize='small')


axs[1].scatter(g_minus_rs, g,marker='.',color='indianred')
axs[1].plot(predicted_x,huber_y,color='seagreen',label=('Huber; ')+r'$\epsilon = $'+str(epsilon_value))
axs[1].invert_yaxis()
axs[1].legend(loc='upper right')

plt.show()