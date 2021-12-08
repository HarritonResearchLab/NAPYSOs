def skyview_plot(image_file,key,center,p_width,deg_width):
    # Import(s)
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import numpy as np
    import pandas as pd
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    # Action

    
    # Get and transform coords
    df = pd.read_csv(key)
    #raw_ras = list(df['RA'])
    #raw_decs = list(df['DEC'])  

    ras = np.array([316.2327596,313.3273727,313.3114793,312.5205471])
    decs = np.array([43.9278512,45.1816813,44.3872614,44.0593017])
    


 
    center_c = SkyCoord(center,unit=(u.deg,u.deg))
    x_cent = float(center_c.ra.deg)
    y_cent = float(center_c.dec.deg)

    #ra_mask = np.logical_and(ras>(x_cent-deg_width/2),ras<(x_cent+deg_width/2))
    #ras = ras[ra_mask]
    #decs = decs[ra_mask]
    #decs = decs[np.logical_and(decs>(y_cent-deg_width/2),decs<(y_cent+deg_width/2))]
    #print(min(ras),max(ras))
    #print(min(decs),max(decs))
    

    # Get parameters to scale coordinates to image

    m = (p_width)/deg_width
    b = (p_width/2)

    left_mask = np.argwhere(ras>x_cent)
    left_ras_prime = (m*(x_cent-ras[left_mask]))+b+134.3
    left_decs_prime = (m*(decs[left_mask]-y_cent))+b+2.7

    right_mask = np.invert(left_mask)
    right_ras_prime = (m*(x_cent-ras[right_mask]))+b-123.6
    right_decs_prime = (m*(decs[right_mask]-y_cent))+b+2.7


    # Plot
    plt.rcParams['font.family'] = 'serif'
    image_data = fits.getdata(image_file,ext=0)
    plt.imshow(image_data,cmap='gray')
    plt.scatter(left_ras_prime,left_decs_prime,color='orange',marker='*',s=5)
    plt.scatter(right_ras_prime,right_decs_prime,color='orange',marker='*',s=5)
    plt.gca().invert_yaxis()
    plt.title('The American and Pelican Nebulae')
    plt.show()

#skyview_plot(image_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/delta/sky_plot/dss2 red ngc7000/dss 2 red 3.2 deg fov.fits',key='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/key.csv',center='314.696+44.33',p_width=1000,deg_width=3.2)

def playing_with_coords(raw_coord):
    #Imports
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    #Action
    c = SkyCoord(raw_coord,unit=(u.deg,u.deg))
    print(c.ra.hourangle)

#playing_with_coords(raw_coord='314.72174296798727+43.9758002436749')

def plot_only(image_file, key):
    #IMports
    from astropy.wcs import WCS
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    # Action

    df = pd.read_csv(key)
    ras = np.array(df['RA'])
    decs = np.array(df['DEC'])

    # Plot
    plt.rcParams['font.family'] = 'serif'
    hdu = fits.open(image_file)[0]
    wcs = WCS(hdu.header)
    plt.subplot(projection=wcs)
    plt.imshow(hdu.data,cmap='gray')
    plt.grid(color='white',ls='solid')
    plt.xlabel('Galactic Longitude')
    plt.ylabel('Galactic Latitude')
    overlay = plt.gca().get_coords_overlay('fk5')
    overlay.grid(color='white',ls='dotted')
    overlay[0].set_axislabel('Right Ascension (J2000)')
    overlay[1].set_axislabel('Declination (J2000)')
    plt.scatter(ras, decs, transform=plt.gca().get_transform('fk5'), s=4, color='orange', edgecolor='black',marker='o')
    plt.show()

plot_only(image_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/delta/sky_plot/best_dss/skyview.fits',key='/home/thaddaeus/FMU/HRL/LAH2.0/data/695/no_clusters/key.csv')
