import numpy as np


def jades_trilegalisl(paths, data, tri_gaia, tri_mag, max_mag, save_file, flag_method, tri_type):
    fieldnum = data.header.fieldnum

    # step 1: compute the area of a JADES image and pull out that many stars
    # from the list at random
    surveyarea = num_pixels * (PIXAR_A2/ 3600**2)


    # step 2: read in the trilegal information

    nfiles = len(V)

    isltot = np.zeros(nfiles)
    islmasked = np.zeros(nfiles)
    islghost = np.zeros(nfiles)

    for jfile in range(0, nfiles):
        # If using trilegal catalog based on UBVRI magnitudes
        if tri_type == 'ubvri':
            mag = V[jfile - 1].mlf

            # step 3: convert from mag to flux
            # Flux calculation
            Fcat = 10 ** (m / (-2.5)) / PIXAR_SR * vega_sirius_Jy


            # step 4: make a mask function

            whpl = V[jfile - 1].V > Max_Mag


        # If using trilegal catalog based on G magnitudes


        # step 5: convert to surface brightness
        magcat = mag[whpl]

        #lIltot = 1e-26 * 1e9 * data.cal.nu * Fcat / (surveyarea * (np.pi / 180) ** 2)
        lIlcat = 1e-26 * 1e9 * data.cal.nu * Fcat[whpl] / (surveyarea * (np.pi / 180) ** 2) # this one is ISL - per radian squared

        #isltot[jfile - 1] = np.sum(lIltot)
        islmasked[jfile - 1] = np.sum(lIlcat)  # this is ISL

    islout['isltotmean'] = np.mean(isltot)
    islout['isltoterr'] = np.std(isltot)

    islout['islmaskedmean'] = np.mean(islmasked)  # this is ISL
    islout['islmaskederr'] = np.std(islmasked)

    data['isl']['tritotmean'] = islout['isltotmean']
    data['isl']['tritoterr'] = islout['isltoterr']
    data['isl']['trimean'] = islout['islmaskedmean']  # this is used as ISL subtraction
    data['isl']['trierr'] = islout['islmaskederr']

    return data



