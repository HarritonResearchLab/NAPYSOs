import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

speedOfLight = 29979200000 #cm/s

#Functions to calculate and plot sound velocities and stellar parameter dependent velocities
#e.g. keplerian and escape velocity
def starParamBased(mass,rOut,starRadius):
    global rStar 
    rStar = starRadius
    def keplerian():
        #mass should be in terms of Solar Masses, rOut should be
        #outer disk radius in AU. 
        global kepVelocity #declaring global variable so it can be called later
        global rOuter 
        global mStar 
        mStar = mass
        rOuter = rOut
        g = 6.7384*(10**-11)
        massStar = mass*1.989*(10**30) #mass of star in kg
        rOutMeters = rOut*1.496*(10**11) #Rout in meters
        kepVelocity = (((massStar*g)/rOutMeters)**(0.5)) #m/s
        kepVelocity = kepVelocity*100 #cm/s
        kepVelocity = round(kepVelocity,0) #correct sig figs
        radiusStar = rStar*69600000000 #radius of the star in cm
        rSunTime = radiusStar/kepVelocity
        timeArray = np.linspace(0, rSunTime,4)
        plt.plot((timeArray*kepVelocity),timeArray,linewidth=5,linestyle='--',alpha=0.5,color='green',label='Keplerian: '+str(kepVelocity)+' cm/s')
        timeArray = np.linspace(rSunTime,maxTime,4)
        plt.plot((timeArray*kepVelocity),timeArray,linewidth=5,linestyle='-',alpha=0.5,color='green')
    def escapeSpeed():
        #radius star should be in terms of stellar radii
        g = 6.7384*(10**-11)
        massStar = mass*1.989*(10**30) #mass of star in kg
        escapeVelocity = ((2*g*massStar)/(starRadius*695500000))**0.5 #m/s
        escapeVelocity = escapeVelocity*100 #cm/s
        escapeVelocity = round(escapeVelocity,0)
        
        radiusStar = rStar*69600000000 #radius of the star in cm
        rSunTime = radiusStar/escapeVelocity
        timeArray = np.linspace(0,rSunTime,10)
        plt.plot((timeArray*escapeVelocity),timeArray,linewidth=5,linestyle='--',alpha=0.5,color='red',label='Escape Velocity; '+str(escapeVelocity)+' cm/s')
        timeArray = np.linspace(rSunTime,maxTime,10)
        plt.plot((timeArray*escapeVelocity),timeArray,linewidth=5,linestyle='-',alpha=0.5,color='red')

    escapeSpeed()
    keplerian()
def soundVelocity(flag,temp,rOut):
    #temp needs to be in Kelvin
    #if flag == '1', temperature is used directly (set rOut 'na'), if flag == '2',
    #then temperature is calculated based on rOut (set temp to 'na')

    if flag == 1:
        if temp > 0:
            speedOfSound = 100*20.05*(temp**0.5) #cm/s
            speedOfSound = round(speedOfSound,0)
            radiusStar = rStar*69600000000 #radius of the star in cm
            rSunTime = radiusStar/speedOfSound
            timeArray = np.linspace(0, rSunTime,4)
            plt.plot((timeArray*speedOfSound),timeArray,linewidth=5,linestyle='--',alpha=0.5,color='blue',label='Sound: '+str(temp)+' K; '+str(speedOfSound)+' cm/s')
            timeArray = np.linspace(rSunTime,maxTime,4)
            plt.plot((timeArray*speedOfSound),timeArray,linewidth=5,linestyle='-',alpha=0.5,color='blue')
        else:
            sys.exit('Error: Temperature must be > 0. Use Kelvin.')
    elif flag == 2:
        tempRouter = rOut
        temp = 300/(rOut**0.75)
        if temp > 0:
            speedOfSound = 100*20.05*(temp**0.5) #cm/s
            speedOfSound=round(speedOfSound,0)
            radiusStar = rStar*69600000000 #radius of the star in cm
            rSunTime = radiusStar/speedOfSound
            timeArray = np.linspace(0, rSunTime,5)
            plt.plot((timeArray*speedOfSound),timeArray,linewidth=5,linestyle='--',alpha=0.5,color='blue',label='Sound: '+str(rOut)+' AU; '+str(round(temp,0))+' K; '+str(speedOfSound)+' cm/s')
            timeArray = np.linspace(rSunTime,maxTime,5)
            plt.plot((timeArray*speedOfSound),timeArray,linewidth=5,linestyle='-',alpha=0.5,color='blue')
        else:
            sys.exit('Error: Temperature must be > 0; rOut must be > 0 as well.')

maxTime = (100*365*30.4167*24*60*60) #number of seconds in 100 years
#timeArray = np.linspace(0,maxTime*2.5,1000) #generates 1000 equally spaced values starting at 0 and ending at max time


#plot
nimbusRoman = {'fontname':'Nimbus Roman'}
plt.style.use('seaborn-darkgrid')

#Dynamic Velocities: 
starParamBased(0.31,290,1.47)
soundVelocity(1,100,'na')
soundVelocity(2,'na',30)
#Constant (user can't change):
radiusStar = rStar*69600000000 #radius of the star in cm
rSunTime = radiusStar/speedOfLight
timeArray = np.linspace(0,rSunTime,4)
plt.plot((timeArray*speedOfLight),timeArray,linewidth=5,linestyle='--',alpha=0.5,color='orange',label='Light')
timeArray = np.linspace(rSunTime,maxTime,4)
plt.plot((timeArray*speedOfLight),timeArray,linewidth=5,alpha=0.5,linestyle='-',color='orange')

#some plot settings
plt.ylim(maxTime,600)
plt.xlim((10*60),(3.08567758*(10**19)))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distance (cm)',**nimbusRoman,fontsize=11)
plt.ylabel('Time (s)',**nimbusRoman,fontsize=11)
plt.yticks(ticks=[(10*60),(60*60),(24*60*60),(30.4368499*24*60*60),(365*30.4368499*24*60*60),(maxTime)],labels=['10 Min','Hour','Day','Month','Year','100 Years'],fontsize=8,fontname='Nimbus Roman')
plt.xticks(ticks=[(rStar*69550000000),(14959787069100),(rOuter*14959787069100),(3085677820000000000)],labels=[r'$R_{star}$',r'$R_{out}$','AU','Parsec'],fontsize=8,fontname='Nimbus Roman')
plt.legend(prop={'family':'Nimbus Roman'})
plt.show()
