def make_plot():
    # Import(s)
    import numpy as np
    import matplotlib.pyplot as plt
    
    # Action
    x =np.linspace(0,15,30)
    y = -0.89*x + np.random.normal(0,1,len(x))
    
   
    
    # Plot
    plt.style.use('seaborn-darkgrid')
    plt.rcParams['font.family'] = 'serif'
    plt.scatter(x,y,color='#408ee0',edgecolor='black',label='Measured Slope')
    plt.xlabel('x')
    plt.ylabel('y')
   
    x1 = plt.gca().get_xlim()[0]
    y1 = plt.gca().get_ylim()[0]
    
    pos_x = x1 - 0.75*np.std(x)
    pos_y = y1 + 0.75*np.std(y)
    plt.scatter(pos_x,pos_y,color='indianred',edgecolor='black',label='Typical Errors')
    
    plt.gca().invert_yaxis()
    
    plt.legend()
    
    plt.show()
    
make_plot()