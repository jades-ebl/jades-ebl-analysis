from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.visualization import PercentileInterval, ImageNormalize
from astropy.visualization import LinearStretch, LogStretch
import os

# Set a flag
flag = 1

# Set the directory containing the FITS files
fits_directory = '/mnt/DataDrive/tsymons/jadez'

# Loop through the FITS files in the directory except catalogue
for filename in os.listdir(fits_directory):
    if filename.endswith('drz.fits'):
        # Construct the full file path
        file_path = os.path.join(fits_directory, filename)

        # Open the FITS file
        with fits.open(file_path) as hdul:



            if flag == 1:
                # Plot the data
                sci_data = hdul[1].data
                plt.figure(figsize=(15, 15))
                #plt.subplot()
                norm = ImageNormalize(sci_data, interval=ZScaleInterval(), stretch=LinearStretch())
                plt.imshow(sci_data, norm=norm, interpolation='none', origin='lower')
                plt.title(filename, fontsize = 20)
                plt.tick_params(axis='both', which='major', labelsize = 18)
                cbar = plt.colorbar()
                cbar.ax.tick_params(labelsize = 18)
                plt.show()