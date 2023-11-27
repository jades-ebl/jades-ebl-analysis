def wcs_getpole(wcsparams):
    import numpy as np
    # define angle constants
    radeg = 180.0 / np.pi
    alpha_0 = wcsparams.crval[0] / radeg
    delta_0 = wcsparams.crval[1] / radeg

    # if theta0 is 90 then this is simple
    if wcsparams.theta0 == 90:
        alpha_p = alpha_0
        delta_p = delta_0
        return alpha_p, delta_p

    # Longpole is the longitude in the native system of the North Pole in the
    # standard system (default = 180 degrees).
    phi_p = wcsparams.lonpole / radeg
    theta_p = wcsparams.latpole / radeg
    sp = np.sin(phi_p)
    cp = np.cos(phi_p)
    sd = np.sin(delta_0)
    cd = np.cos(delta_0)
    tand = np.tan(delta_0)

    # If origin is set then crval gives the coordinates of the origin in the
    # native system.   This must be converted (using Eq. 7 in
    # Greisen & Calabretta with theta0 = 0) to give the coordinates
    # of the North pole (alpha_p, delta_p)
    if wcsparams.theta == 0:
        if delta_0 == 0 and np.abs(wcsparams.lonpole) == 90.0:
            delta_p = theta_p
        else:
            delta_p = np.arccos(sd / cp)

        if wcsparams.latpole != 90 and np.abs(theta_p + delta_p) < np.abs(theta_p - delta_p):
            delta_p = -delta_p

        if wcsparams.lonpole == 1.8e2 or cd == 0.0:
            alpha_p = alpha_0
        else:
            alpha_p = alpha_0 - np.arctan2(sp / cd, -np.tan(delta_p) * tand)
    else:
        ctheta = np.cos(wcsparams.theta0 / radeg)
        stheta = np.sin(wcsparams.theta0 / radeg)
        term1 = np.arctan2(stheta, ctheta * cp)
        term2 = np.arccos(sd / np.sqrt(1.0 - ctheta ** 2 * sp ** 2))

        if term2 == 0.0:
            delta_p = term1
        else:
            delta_p1 = np.abs((term1 + term2) * radeg)
            delta_p2 = np.abs((term1 - term2) * radeg)

            if delta_p1 > 90 and delta_p2 > 90:
                print('wcs_getpole: No valid solution!')

            if delta_p1 <= 90 and delta_p2 > 90:
                delta_p = term1 + term2

            if delta_p1 > 90 and delta_p2 <= 90:
                delta_p = term1 - term2

            if delta_p1 <= 90 and delta_p2 <= 90:
                # there are two valid solutions
                delta_p1 = (term1 + term2) * radeg
                delta_p2 = (term1 - term2) * radeg

                if np.abs(wcsparams.latpole - delta_p1) < np.abs(wcsparams.latpole - delta_p2):
                    delta_p = term1 + term2
                else:
                    delta_p = term1 - term2

            if cd == 0.0:
                alpha_p = alpha_0
            else:
                sdelt = np.sin(delta_p)

                if sdelt == 1:
                    alpha_p = alpha_0 - phi_p - np.pi
                elif sdelt == -1:
                    alpha_p = alpha_0 - phi_p
                else:
                    alpha_p = alpha_0 - np.arctan2((stheta - np.sin(delta_p) * sd) / (np.cos(delta_p) * cd),
                                                   sp * ctheta / cd)

    return alpha_p, delta_p
