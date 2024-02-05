import numpy as np
from scipy.io import loadmat

def nh_trilegalisl(paths, data, tri_gaia, tri_mag, max_mag, save_file, flag_method, tri_type):
    fieldnum = data.header.fieldnum

    # step 1: compute the area of a LORRI(JADES) image and pull out that many stars
    # from the list at random
    surveyarea = (data.astrom.imagew
                         * data.astrom.imageh
                         * data.cal.pixsize_arcsec ** 2) / 3600 ** 2
    # num of pix in image * degree/pix = total deg^2


    # step 2: read in the trilegal information
    # determine which files are being examined
    # path??
    if paths.datadir == '/data/symons/NH_old_data/mat/ghosts/':
        ghost = 1
        old = 0
        new = 0
    elif paths.datadir == '/data/symons/NH_old_data/mat/good/':
        old = 1
        ghost = 0
        new = 0
    elif paths.datadir == '/data/symons/nh_data/mat/':
        new = 1
        old = 0
        ghost = 0
    elif paths.datadir == '/data/symons/nh_data_lauer/mat/':
        new = 0
        old = 0
        ghost = 0
        lauer = 1
    elif paths.datadir == '/data/symons/nh_data_new/mat/':
        new = 0
        old = 0
        ghost = 0
        lauer = 0
        newest = 1

    # load the appropriate trilegal catalog files
    # path??
    if new == 1:
        filename = f'.mat'
        trilegal_data = loadmat(filename)
    elif old == 1:
        filename = f'.mat'
        trilegal_data = loadmat(filename)
    elif lauer == 1:
        filename = f'.mat'
        trilegal_data = loadmat(filename)
    elif newest == 1:
        filename = f'.mat'
        trilegal_data = loadmat(filename)

    nfiles = len(V)

    isltot = np.zeros(nfiles)
    islmasked = np.zeros(nfiles)
    islghost = np.zeros(nfiles)

    for jfile in range(0, nfiles):
        # If using trilegal catalog based on UBVRI magnitudes
        if tri_type == 'ubvri':
            mag = V[jfile - 1].mlf

            # step 3: convert from LORRI-band mag to flux
            if flag_method == 'old_corr' or flag_method == 'new':
                Fcat = data.cal.vzero * 10 ** (-V[jfile - 1].mlf / 2.5)
            elif flag_method == 'old':
                Fcat = 3055 * 10 ** (-V[jfile - 1].mlf / 2.5)

            # step 4: make a mask function
            if tri_gaia == 0:
                if flag_method == 'old_corr' or flag_method == 'new':
                    whpl = V[jfile - 1].V > max_mag
                    # produces a boolean array
                elif flag_method == 'old':
                    whpl = V[jfile - 1].V > data.mask.maxmag
            elif tri_gaia == 1
                whpl = V[jfile - 1].V > tri_mag

        # If using trilegal catalog based on G magnitudes
        elif tri_type == 'gaia':
            mag = V[jfile - 1].G

            # step 3: convert from LORRI-band mag to flux
            if flag_method == 'old_corr' or flag_method == 'new':
                Fcat = data.cal.vzero * 10 ** (-V[jfile - 1].G / 2.5)
            elif flag_method == 'old':
                Fcat = 3055 * 10 ** (-V[jfile - 1].G / 2.5)

            # step 4: make a mask function
            if tri_gaia == 0:
                if flag_method == 'old_corr' or flag_method == 'new':
                    whpl = V[jfile - 1].G > max_mag
                elif flag_method == 'old':
                    whpl = V[jfile - 1].G > data.mask.maxmag
            elif tri_gaia == 1:
                whpl = V[jfile - 1].G > tri_mag

        # step 5: convert to surface brightness
        magcat = mag[whpl]

        lIltot = 1e-26 * 1e9 * data.cal.nu * Fcat / (surveyarea * (np.pi / 180) ** 2)
        lIlcat = 1e-26 * 1e9 * data.cal.nu * Fcat[whpl] / (surveyarea * (np.pi / 180) ** 2) # this one is ISL - per radian squared

        isltot[jfile - 1] = np.sum(lIltot)
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



