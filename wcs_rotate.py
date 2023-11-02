import numpy as np

def wcs_rotate(wcsparams):
    # define angle constants
    radeg = 1.8e2 / np.pi

    # check that the input fields make sense
    if 'phi' in wcsparams or 'theta' in wcsparams:
        if 'phi' not in wcsparams or 'theta' not in wcsparams:
            print("wcs_rotate error: 'theta' and 'phi' are not set correctly! Aborting.\n")
            return
        reverseme = 1
    else:
        if 'longitude' not in wcsparams or 'latitude' not in wcsparams:
            print("wcs_rotate error: 'longitude' and 'latitude' are not set correctly! Aborting.\n")
            return
        reverseme = 0

    if 'theta0' not in wcsparams:
        print("wcs_rotate error: 'theta0' not set! Aborting.\n")

    # Longpole is the longitude in the native system of the North Pole
    # in the standard system (default = 180 degrees).
    if 'lonpole' not in wcsparams:
        wcsparams['lonpole'] = 1.8e2

    phi_p = wcsparams['lonpole'] / radeg
    sp = np.sin(phi_p)
    cp = np.cos(phi_p)


    # If Theta0 = 90 then CRVAL gives the coordinates of the origin in the
    # native system.   This must be converted (using Eq. 7 in Greisen & Calabretta
    # with theta0 = 0) to give the coordinates of the North pole.
    alpha_p, delta_p = wcs_getpole(wcsparams)

    # compute useful quantities relating to reference angles
    sa = np.sin(alpha_p)
    ca = np.cos(alpha_p)
    sd = np.sin(delta_p)
    cd = np.cos(delta_p)

    # calculate rotation matrix
    r = np.array([[-sa * sp - ca * cp * sd, ca * sp - sa * cp * sd, cp * cd],
                  [sa * cp - ca * sp * sd, -ca * cp - sa * sp * sd, sp * cd],
                  [ca * cd, sa * cd, sd]])

    # solve the set of equations for each datum point
    if reverseme:
        latitude = wcsparams['phi']
        longitude = wcsparams['theta']
        g = np.where(np.isfinite(wcsparams['phi']) & np.isfinite(wcsparams['theta']))

        if len(g) == 0:
            # This would be an error
            print("wcs_rotate: No finite input theta, phi given!")
            return
        else:
            phi1 = wcsparams['phi'][g] / radeg
            theta1 = wcsparams['theta'][g] / radeg
            r = r.T
    else:
        phi = wcsparams['longitude']
        phi1 = wcsparams['longitude'] / radeg
        theta1 = wcsparams['latitude'] / radeg

    # Define the right-hand side of the equations
    l = np.cos(theta1) * np.cos(phi1)
    m = np.cos(theta1) * np.sin(phi1)
    n = np.sin(theta1)

    # find solution to the system of equations and put it in b
    # can't use matrix notation in case l,m,n are vectors
    b0 = np.sum(r[0] * [l, m, n], axis=0)
    b1 = np.sum(r[1] * [l, m, n], axis=0)
    b2 = np.sum(r[2] * [l, m, n], axis=0)

    if reverseme:
        wcsparams['latitude'][g] = np.arcsin(b2) * radeg
        wcsparams['longitude'][g] = np.arctan2(b1, b0) * radeg
    else:
        wcsparams['theta'] = np.arcsin(b2) * radeg
        wcsparams['phi'] = np.arctan2(b1, b0) * radeg

    return wcsparams

