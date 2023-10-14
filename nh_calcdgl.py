def nh_calcdgl(data, paths, flag_method):
    dglparams = nh_get_dgl_params()
    # calls another function, "nh_get_dgl_params", and assigns the output to "dglparams"
    # to be defined in another script

    import os
    from astropy.io import fits
    import scipy.io as sp
    import numpy as np

    want = 'iris'

    # compares the variable "want" to the string "planck"
    # returns 0 if false, 1 if true
    if want == 'planck':

        # Checks if Planck map exists, if not calls script that makes them
        planck_filename = os.path.join(paths.planckdir, f'planck_{data.header.timestamp}_fx.fits')

        # If at least one file is not found, calls script to make them all
        if not os.path.isfile(planck_filename):
            print('No Planck file found, retrieving new Planck file.')
            pydir = '/home/symons/nh_ebl_pipeline/py/Planck-Cirrus_Estimation-master/'
            pyfile = 'get_planck.py'
            imagefile = os.path.join(paths.imagedir, f"regist_{data.header.rawfile}")

            # Writes the imagefile path in a text file to where python script is
            with open(os.path.join(pydir, 'imagefile.txt'), 'w') as file:
                file.write(imagefile)

                os.system(f'python {os.path.join(pydir, pyfile)}')

                # Moves the Planck files to the data directory
                if not os.path.isdir(paths.planckdir):
                    os.mkdir(paths.planckdir)

                os.rename(os.path.join(pydir, f'planck_{data.header.timestamp}_fx.fits'),
                          os.path.join(paths.planckdir, f"planck_{data.header.timestamp}_fx.fits"))

        # Constructs FITS file path
        planck_fits_file_path = planck_filename

        # Reads the FITS file
        with fits.open(planck_fits_file_path) as hdul:
            irismap = hdul[0].data

        irisim = irismap

        type = 'planck'

    # Check if planck maps are already made, if not call python script that makes
    elif want == 'planck_mc':

        # If at least one of the files isn't there, call python script to make them all
        planck_mc_filename = os.path.join(paths.panckmcdir, f'planck_{data.header.timestamp}_errmean.mat')

        if not os.path.isfile(planck_mc_filename):
            print('No Planck MC file found, retrieving new Planck MC file.')
            pydir = '/home/symons/nh_ebl_pipeline/py/Planck_Cirrus_Estimation-master/'
            pyfile = 'get_planck_mc.py'
            imagefile = os.path.join(paths.imagedir, f'regist_{data.header.rawfile}')

            # Write the imagefile path in a text file to where the python script is
            with open(os.path.join(pydir, 'imagefile.txt'), 'w') as file:
                file.write(imagefile)

            # Call the Python script
            os.system(f'python {os.path.join(pydir, pyfile)}')

            # Create the necessary directory if it doesn't exist
            if not os.path.exists(paths.planckmcdir):
                os.mkdir(paths.planckmcdir)

            # Move the Planck MC files from the Python directory to the specified directory
            os.rename(
                os.path.join(pydir, f'planck_{data.header.timestamp}_errmean.mat'),
                os.path.join(paths.planckmcdir, f'planck_{data.header.timestamp}_errmean.mat'))

        # Load the MAT-file
        mat_file_path = os.path.join(paths.planckmcdir, f'planck_{data.header.timestamp}_errmean.mat')
        loaded_data = sp.loadmat(mat_file_path)

        err_means = loaded_data['err_means']

        # Tau = 1, Temperature = 2, Beta = 3 | Opacity = 1, Temperature = 2, Spectral-Index = 3
        # Defines constants for indices
        TAU_INDEX = 0
        TEMP_INDEX = 1
        BETA_INDEX = 2

        # Extract data and perform calculations: tau
        planck_mc_100m_tau = err_means[TAU_INDEX, :]
        planck_mc_mean_tau = np.mean(err_means[TAU_INDEX, :])
        planck_mc_sem_tau = np.std(err_means[TAU_INDEX, :]) / np.sqrt(len(planck_mc_100m_tau))

        # Temp
        planck_mc_100m_temp = err_means[TEMP_INDEX, :]
        planck_mc_mean_temp = np.mean(err_means[TEMP_INDEX, :])
        planck_mc_sem_temp = np.std(err_means[TEMP_INDEX, :]) / np.sqrt(len(planck_mc_100m_temp))

        # Beta
        planck_mc_100m_beta = err_means[BETA_INDEX, :]
        planck_mc_mean_beta = np.mean(err_means[BETA_INDEX, :])
        planck_mc_sem_beta = np.std(err_means[BETA_INDEX, :]) / np.sqrt(len(planck_mc_100m_beta))

        type = 'planck_mc'

    # Rewrite code without subprocess
    elif want == 'iris':

        # Check if iris maps are already made, if not call python script that makes them
        # Define FITS file path for IRIS
        iris_fits_file_path = os.path.join(paths.irisdir, f'iris_{data.header.timestamp}_fx.fits')

        # If at least one of the files isn't there, call python script that makes them all
        if not os.path.isfile(iris_fits_file_path):
            print('No IRIS file found, retrieving new IRIS file.')

            pydir = '/home/symons/nh_ebl_pipeline/py/IRISpy-master/'  # Directory where the python script is
            pyfile = 'nh_iris_map.py'  # Name of python file to run
            imagefile = os.path.join(paths.imagedir,
                                     f'regist_{data.header.rawfile}')  # FITS file to send to python script

            # Writes the imagefile path in a text file to where the python script is
            with open(os.path.join(pydir, 'imagefile.txt'), 'w') as file:
                file.write(imagefile)

            # Call the Python script
            os.system(f'python {os.path.join(pydir, pyfile)}')

            # Creates the IRIS directory if it doesn't exist
            if not os.path.isdir(paths.irisdir):
                os.makedirs(paths.irisdir)

            # Gets IRIS files and moves them to their data directory
            iris_source_file = os.path.join(pydir, f'iris_{data.header.timestamp}_fx.fits')
            iris_destination_file = os.path.join(paths.irisdir, f'iris_{data.header.timestamp}_fx.fits')

            # Moves from pydir to paths.irisdir
            if os.path.isfile(iris_source_file):
                os.rename(iris_source_file, iris_destination_file)

        # Load in IRIS map
        irismap = fits.getdata(iris_fits_file_path)  # Reads FITS file

        irisim = irismap

        type = 'iris'

    elif want == 'iris_sfd':

        # Checks if IRIS SFD maps are already made, if not call python script that makes them
        # Defines FITS file path for IRIS SFD
        iris_sfd_fits_file_path = os.path.join(paths.irissfddir, f'iris_sfd_{data.header.timestamp}_fx.fits')

        if not os.path.isfile(iris_sfd_fits_file_path):
            print('No IRIS SFD file found, retrieving new IRIS SDF file.')

            pydir = '/home/symons/nh_ebl_pipeline/py/IRIS_SFD_Estimation/'  # Directory where python script is
            pyfile = 'get_iris_sfd.py'  # Name of python file to run
            imagefile = os.path.join(paths.imagedir,
                                     f'regist_{data.header.rawfile}')  # FITS file to send the python script

            # Write the imagefile path in a text file where the python script is
            with open(os.path.join(pydir, 'imagefile.txt'), 'w') as file:
                file.write(imagefile)

            # Call the Python script
            os.system(f'python {os.path.join(pydir, pyfile)}')

            # Create the IRIS SFD directory if it doesn't exist
            if not os.path.isdir(paths.irissfddir):
                os.makedirs(paths.irissfddir)

            # Gets the generated IRIS SFD file and moves it to its data directory
            iris_sfd_source_file = os.path.join(pydir, f'iris_sfd_{data.header.timestamp}_fx.fits')
            iris_sfd_destination_file = os.path.join(paths.irissfddir, f'iris_sfd_{data.header.timestamp}_fx.fits')

            # Moves from pydir to paths.irissfddir
            if os.path.isfile(iris_sfd_source_file):
                os.rename(iris_sfd_source_file, iris_sfd_destination_file)

        # Load in planck maps
        irismap = fits.getdata(iris_sfd_fits_file_path)  # Reads FITS file

        irisim = irismap

        type = 'iris_sfd'

    elif want == 'nh':

        nh_fits_file_path = os.path.join(paths.nhdir, f'nh_{data.header.timestamp}_fx.fits')

        if not os.path.isfile(nh_fits_file_path):
            print('No NH file found, retrieving new NH file.')

            pydir = '/home/symons/nh_ebl_pipeline/py/NH_Estimation/'  # Directory where python script is
            pyfile = 'get_nh.py'  # Name of python file to run
            imagefile = os.path.join(paths.imagedir,
                                     f'regist_{data.header.rawfile}')  # FITS file to send the python script

            # Write the imagefile path in a text file where the python script is
            with open(os.path.join(pydir, 'imagefile.txt'), 'w') as file:
                file.write(imagefile)

            # Call the Python script
            os.system(f'python {os.path.join(pydir, pyfile)}')

            # Create the NH directory if it doesn't exist
            if not os.path.isdir(paths.nhdir):
                os.makedirs(paths.nhdir)

            # Gets the generated NH file and moves it to its data directory
            nh_source_file = os.path.join(pydir, f'nh_{data.header.timestamp}_fx.fits')
            nh_destination_file = os.path.join(paths.nhdir, f'nh_{data.header.timestamp}_fx.fits')

            # Moves from pydir to paths.nhdir
            if os.path.isfile(nh_source_file):
                os.rename(nh_source_file, nh_destination_file)

        # Loads in NH maps
        irismap = fits.getdata(nh_fits_file_path)

        irisim = irismap

        type = 'nh'

    if type in ['iris', 'iris_sfd', 'nh', 'planck']:
        if type == 'iris':

            # Calls external function "nh_dgl_manmask" to create manual masks (man masks)
            data = nh_dgl_manmask(data, paths)

            # Load in man masks
            filein = f'{paths.dglmandir}{data.header.timestamp}.mat'

            # If man mask files exists, create the mask, otherwise no man mask
            if os.path.exists(filein):
                manmask = np.load(filein)
            else:
                manmask = np.zeros(data.data.shape)  # Won't be used, empty masks are saved for all

                # Calculate mean and std of iris map
                ohm_mean = np.mean(irisim)
                ohm_std = np.std(irisim)

            # If not looking at Neutral Hydrogen, calculate DGL, otherwise skip
            if type == 'nh':
                if type == 'iris':

                    # nu I_nu is (iris map - 0.8 MJy/sr) * 1e-20 * 3e8 / 100e-6 * 1e9
                    onehundo = np.mean(irisim[~manmask] - 0.48)  # mean 100 micron emission including CIB subtraction
                    # (real value 0.48, testing 0.24)

                    nuinu = (irisim[~manmask] - 0.48) * dglparams.norm  # Old value was 0.8 MJy/sr, now using 0.48
                    # MJy/sr from 14.4 +/- 6.3 nW in https://articles.adsabs.harvard.edu/pdf/2011ASPC..446..309B

                elif type == 'iris_sfd':

                    # nu I_nu is (iris map - 0.8 MJy/sr) * 1e-20 * 3e8 / 100e-6 * 1e9
                    onehundo = np.mean(irisim)  # mean 100 micron emission including CIB subtraction
                    nuinu = irisim * dglparams.norm  # IRIS/SFD doesn't need CIB sub if zero-referenced to NHI

                elif type == 'planck':

                    # nu I_nu is (iris map) * 1e-20 * 3e8/100e-6 * 1e9 - planck map already
                    onehundo = np.mean(irisim)  # mean 100 micron emission including CIB subtraction
                    nuinu = irisim * dglparams.norm

                # b_lambda (I_nu_opt/I_nu_100um) is estimated by fitting mean ZDA04 model to many measurements
                # of b_lambda
                # c_lambda (lambda*I_lambda_opt/nu*I_nu_100um)s is 10^-6*(100um/0.655um)*b_lambda
                # cbar is 0.491 (cbar_655nm calculated by integrating cbar over lorri bandpass)
                cbar = dglparams.cbar[0]

                # A = 1/0.567 is normalizing factor d0 (from normalizing d(b) at b=25)
                # g = 0.61 is bandpass-weighted mean of observation-constrained model for
                # high-lat diffuse dust component of DGL
                # dl is d(b), a function that accounts for change in c_lambda due to b
                dl = dglparams.A * (1 - 1.1 * dglparams.g[0] * np.sqrt(np.sin(np.abs(data.coords.galactic[1]) *
                                    np.pi / 180)))

                # dglim is estimate of DGL using average of 100um intensity, c_lambda, and d(b)
                dglim = nuinu * cbar * dl

                # propagated error here
                if type == 'iris':
                    dglerr_nuinu = (cbar * dl)**2 * 0.21 * 2  # New CIB error is 0.21 MJy/sr
                    dglerr_cbar = (np.mean(nuinu[~manmask]) * dl)**2 * dglparams.cbar[1]**2
                    dglerr_dl = ((np.mean(nuinu[~manmask])) * cbar * dglparams.A * 1.1 *
                                 np.sqrt(np.sin(np.abs(data.coords.galactic[1]) * np.pi / 180))**2 * dglparams.g[1]**2)

                    dgl_err = np.sqrt(dglerr_nuinu + dglerr_cbar + dglerr_dl)

                else:
                    dglerr_nuinu = (cbar * dl)**2 * 0**2  # No error from CIB if not subtracting
                    dglerr_cbar = (np.mean(nuinu) * dl)**2 * dglparams.cbar[1]**2
                    dglerr_dl = (np.mean(nuinu) * cbar * dglparams.A * 1.1 *
                                 np.sqrt(np.sin(np.abs(data.coords.galactic[1]) * np.pi / 180))**2 * dglparams.g[1]**2)

                    dgl_err = np.sqrt(dglerr_nuinu + dglerr_cbar + dglerr_dl)

                # mean and std of dgl estimate (at each pixel?)
                if type == 'iris':
                    dgl_mean = np.mean(dglim[~manmask])
                    dgl_std = np.std(dglim[~manmask])

                else:
                    dgl_mean = np.mean(dglim)
                    dgl_std = np.std(dglim)

                if type == 'iris':
                    data.dgl.dglmean_iris = dgl_mean
                    data.dgl.dglerr_iris = dgl_err
                    data.dgl.dglstd_iris = dgl_std
                    data.dgl.ohmmean_iris = ohm_mean
                    data.dgl.ohmstd_iris = ohm_std
                    data.dgl.ohmim_iris = irisim * ~manmask
                    data.dgl.onehundomean_iris = onehundo
                    data.dgl.dglim_iris = dglim
                    data.dgl.conv_iris = dglparams.cbar[0] * dl
                    data.dgl.convp_iris = dglparams.norm * dglparams.cbar[0] * dl
                elif type == 'iris_sfd':
                    data.dgl.dglmean_iris_sfd = dgl_mean
                    data.dgl.dglerr_iris_sfd = dgl_err
                    data.dgl.dglstd_iris_sfd = dgl_std
                    data.dgl.ohmmean_iris_sfd = ohm_mean
                    data.dgl.ohmstd_iris_sfd = ohm_std
                    data.dgl.ohmim_iris_sfd = irisim
                    data.dgl.onehundomean_iris_sfd = onehundo
                    data.dgl.dglim_iris_sfd = dglim
                    data.dgl.conv_iris_sfd = dglparams.cbar[0] * dl
                    data.dgl.convp_iris_sfd = dglparams.norm * dglparams.cbar[0] * dl
                elif type == 'planck':
                    data.dgl.dglmean_planck = dgl_mean
                    data.dgl.dglerr_planck = dgl_err
                    data.dgl.dglstd_planck = dgl_std
                    data.dgl.ohmmean_planck = ohm_mean
                    data.dgl.ohmstd_planck = ohm_std
                    data.dgl.ohmim_planck = irisim
                    data.dgl.onehundomean_planck = onehundo
                    data.dgl.dglim_planck = dglim
                    data.dgl.conv_planck = dglparams.cbar[0] * dl
                    data.dgl.convp_planck = dglparams.norm * dglparams.cbar[0] * dl
                elif type == 'nh':
                    data.dgl.ohmmean_nh = ohm_mean
                    data.dgl.ohmstd_nh = ohm_std
                    data.dgl.ohmim_nh = irisim
                elif type == 'planck_mc':
                    data.dgl.ohm_planck_mc_tau = planck_mc_100m_tau
                    data.dgl.mean_planck_mc_tau = planck_mc_mean_tau
                    data.dgl.sem_planck_mc_tau = planck_mc_sem_tau
                    data.dgl.ohm_planck_mc_temp = planck_mc_100m_temp
                    data.dgl.mean_planck_mc_temp = planck_mc_mean_temp
                    data.dgl.sem_planck_mc_temp = planck_mc_sem_temp
                    data.dgl.ohm_planck_mc_beta = planck_mc_100m_beta
                    data.dgl.mean_planck_mc_beta = planck_mc_mean_beta
                    data.dgl.sem_planck_mc_beta = planck_mc_sem_beta
