""" JADES_ISL.py

    Calculates ISL contribution from psf wings

    5/3/24
"""

import numpy as np
from RegridderShiftNumba import RegridderShift_Numba as RegridderShift
import scipy

def input_mapper(ncat,im_shape,ypix,xpix,flux,fluxerr,errflag_psf):
    """input_mapper(ncat,im_shape,ypix,xpix,flux,fluxerr,errflag_psf)

    Args:
        ncat ('int'): length of catalog
        im_shape ('tuple'): shape of data
        ypix ('numpy.ndarray'): • Y - pix (double), y-centroid in detection image, windowed
        xpix ('numpy.ndarray'): • X - pix (double), x-centroid in detection image, windowed
        flux ('numpy.ndarray'): flux column from catalog
        fluxerr ('numpy.ndarray'): corresponded flux error from catalog
        errflag_psf (0 or 1): 0- using flux from catalog, 1- error calculation 

    Returns:
        sourceMap('numpy.ndarray'): sources placed via subpixel shifting into simulated map of zeros
    """

    #seed random numbers used for error calculation
    np.random.seed(1)

    #Create empty lists to append coordinates and flux
    sourceCoordRound_accumulator = []
    sourceCoordDelta_accumulator = []
    sourceCoordFlux_accumulator = []

    # Reference source is an ideal point source with 2 pixels on either side to allow for source shifting beyond integer locations
    sourceRef = np.zeros((1, 1))
    sourceSize = np.int64(sourceRef.shape[0])
    sourceRef[sourceSize // 2, sourceSize // 2] = 1.0

    # Create simulated map (zeros)
    sourceMap = np.zeros(im_shape)

    #In the range of all the sources in catalog, append coordinates and flux    
    for i in range(ncat):
        
        #Center of source in pixels
        y,x= ypix[i], xpix[i]
        sourceCoord = [x,y]
        #Flux of source
        if errflag_psf ==1:
            sourceCoordFlux= flux[i] + (fluxerr[i] * np.random.randn(1)[0])
            
        elif errflag_psf==0:
            sourceCoordFlux= flux[i]
        
        sourceCoordRound = np.int64(np.round(sourceCoord))  # Make a consistent set of rounded coordinates
        sourceCoordDelta = sourceCoord - sourceCoordRound  # Get the decimal difference 
        
        #append values
        sourceCoordRound_accumulator.append(sourceCoordRound)
        sourceCoordDelta_accumulator.append(sourceCoordDelta)
        sourceCoordFlux_accumulator.append(sourceCoordFlux)

    # Place the sources via subpixel shifting
    sourceMap = sourceMapper(ncat,sourceMap, sourceCoordRound_accumulator, sourceCoordDelta_accumulator, sourceCoordFlux_accumulator, sourceRef, sourceSize)
    
    return sourceMap

def sourceMapper(ncat,sourceMap,sourceCoordRound,sourceCoordDelta,sourceCoordFlux,sourceRef,sourceSize):
    """sourceMapper(ncat,sourceMap,sourceCoordRound,sourceCoordDelta,sourceCoordFlux,sourceRef,sourceSize)

    Args:
        ncat ('int'): length of catalog
        sourceMap('numpy.ndarray'): simulated map of zeros
        sourceCoordRound ('list'): rounded coordinates (x,y)
        sourceCoordDelta ('list'): decimal difference between coordinate of source and rounded coordinate
        sourceCoordFlux ('list'): flux of sources
        sourceRef ('numpy.ndarray'): reference source is an ideal point source ( np.zeros((1, 1)) )
        sourceSize ('int64'): size of ideal point source (1)

    Returns:
        sourceMap('numpy.ndarray'): sources placed via subpixel shifting into simulated map of zeros
    """

    # Place the sources via subpixel shifting
    for i in range(ncat):
        # Source size here is the unshifted one (so does not take into account extra row/column added by shifting - that account happens in the if statement set)
        xp = [sourceCoordRound[i][ 0] - np.int64(np.floor(sourceSize / 2)),
              sourceCoordRound[i][0] + np.int64(np.ceil(sourceSize / 2))]
        yp = [sourceCoordRound[i][1] - np.int64(np.floor(sourceSize / 2)),
              sourceCoordRound[i][1] + np.int64(np.ceil(sourceSize / 2))]
        
        # Scale and shift the source
        sourceCurrent = RegridderShift(sourceRef * sourceCoordFlux[i],[sourceCoordDelta[i][ 0], sourceCoordDelta[i][ 1]])

        # Adjust due to shifting by updating the edges
        if (sourceCoordDelta[i][ 0] > 0):
            xp[1] += np.int64(np.ceil(sourceCoordDelta[i][ 0]))
        elif (sourceCoordDelta[i][ 0] < 0):
            xp[0] += np.int64(np.floor(sourceCoordDelta[i][ 0]))
        
        if (sourceCoordDelta[i][1] > 0):
            yp[1] += np.int64(np.ceil(sourceCoordDelta[i][ 1]))
        elif (sourceCoordDelta[i][ 1] < 0):
            yp[0] += np.int64(np.floor(sourceCoordDelta[i][ 1]))
        
        sourceMap[yp[0]:yp[1], xp[0]:xp[1]] += sourceCurrent  # Write in the source

    return sourceMap

def ISL(psf,sourceMap,wavelength,omega_pix,mask):
    """ISL(psf,sourceMap,wavelength,omega_pix,mask):
        #calculates surface brightness of nonmasked psf wings

    Args:
        psf ('numpy.ndarray'): PSF of filter
        sourceMap ('numpy.ndarray'): sources placed via subpixel shifting into simulated map of zeros
        wavelength ('float'): Pivot Wavelength in micrometers
        omega_pix ('float'): Pixel Area in Steradians
        mask ('numpy.ndarray'): array of masked values and 1s

    Returns:
        sourceMap_convolved('numpy.ndarray'): convolved array of noiseless simulation
        ISL_mean('float'): surface brightness of nonmasked psf wings
    """
    #replace non-finite values (nans and infs) with 0
    sourceMap=np.where(np.isfinite(sourceMap)==False,0,sourceMap)
    #convolve sourceMap
    sourceMap_convolved= scipy.signal.fftconvolve(sourceMap,psf,mode='same')
    
    #convert to surface brightness in nW m^-2/str
    surface_brightness= ((sourceMap_convolved* 10**-9)* ((3*10**8)/(wavelength*10**(-6)))* 10**(-26) * 10**9)/omega_pix

    #apply sources mask of 0s and 1s to convolved array
    masked_sim=surface_brightness*mask
    #avoid forced zeros by replacing with nan
    masked_sim=np.where(masked_sim==0,np.nan,masked_sim)

    #mean
    ISL_mean = np.nanmean(masked_sim)


    return sourceMap_convolved, ISL_mean

def ISL_convolve(psf,sourceMap):
    """ISL_convolve(psf,sourceMap)
        #convolves array of noiseless simulation

    Args:
        psf ('numpy.ndarray'): PSF of filter
        sourceMap ('numpy.ndarray'): sources placed via subpixel shifting into simulated map of zeros

    Returns:
        sourceMap_convolved('numpy.ndarray'): convolved array of noiseless simulation
    """
    
    #convolve sourceMap
    sourceMap_convolved= scipy.signal.fftconvolve(sourceMap,psf,mode='same')

    return sourceMap_convolved

