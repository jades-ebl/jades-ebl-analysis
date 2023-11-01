import numpy as np
def wcsxy2sph(longitude, latitude, wcsparams):
    # define angle constants
    radeg = 180 / np.pi

    # find the number of elements in each of the data arrays
    n_long = len(longitude)
    n_lat = len(latitude)
    
    # error checking
    if n_long != n_lat:
        print("wcssph2xy.py: longitude and latitude inputs must have the same number of elements.\n")
        
    if 'lonpole' not in wcsparams:
        print("wcssph2xy.py: wcsparams structure field lonpole missing!\n")
        return
    if 'latpole' not in wcsparams:
        print("wcssph2xy.py: wcsparams structure field latpole missing!\n")
        return
    if 'crval' not in wcsparams:
        print("wcssph2xy.py: wcsparams structure field crval missing!\n")
        return
    else:
    if len(wcsparams['crval']) != 2:
        print("wcssph2xy.py: wcsparams structure field crval must have 2 elements!\n")
        return
    
    # Convert all longitude values into the range -180 to 180 so that equations work properly.
    lng = longitude.copy()
    lat = latitude.copy()

    temp = numpy.where(lng > 180.0)
    Ntemp = len(temp)

    if Ntemp > 0:
        lng[temp] -= 360.0
    
    
    offset = 1.0e-7
    bad = numpy.where(numpy.abs(lat - 90.0) < offset * radeg)

    if len(bad) > 0:
        lat[bad] = 90.0 - north_offset * radeg
        
    bad = np.where(np.abs(lat + 90.0) < offset * radeg)
    if len(bad) > 0:
        lat[bad] = south_offset * radeg - 90.0
        
    # append to wcsparams for passing to rotation function
    wcsparams['longitude'] = lng
    wcsparams['latitude'] = lat
    wcsparams['theta0'] = 90.0
    
    # convert from standard coordinate system to "native" coordinate system
    wcsparams = wcs_rotate(wcsparams)
    
    # strip out the parts of interest
    phi = wcsparams['phi'] / radeg
    theta = wcsparams['theta'] / radeg
    
    # This inversion is only for the TAN projection - more can be added
    # but that requires this get changed.
    xpix = np.nan * np.zeros(theta.shape)
    ypix = xpix
    g = np.where(theta > 0)
    Ng = len(g)
    if Ng > 0:
        r_theta = radeg / np.tan(theta[g])
        xpix[g] = r_theta * np.sin(phi[g])
        ypix[g] = -r_theta * np.cos(phi[g])
    
    
    return xpix, ypix
