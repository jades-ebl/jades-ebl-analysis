
""" MAIN_JADES_driver.py
    
    Use MAIN_JADES_config file in order to set variables
    Reads in JADES image data for filter, catalog file for a version of data, Nircam PSF, and Trilegal Pickle File
    Calculates EBL

    9/19/24
"""

#Imports
from astropy.io import fits
import os
from astropy.table import Table
import astropy.wcs as wcs
from astropy.io import fits
import numpy as np
import pickle
#plotting
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval,LinearStretch,ImageNormalize
#Import configuration variables used to read in image data, catalog data, masking, and histogram
from MAIN_JADES_config import *
#import Masking Function
from JADES_Mask import Final_Masking_new
#Import PSF ISL Function
from JADES_ISL import input_mapper,sourceMapper,ISL
#Import function to calculate magnitude zeropoint
from JADES_calc_zeropoint import calc_zeropoint
#Import faint ISL Function
from JADES_trilegalisl import jades_trilegalisl, jades_trilegalisl_chunk


#Read in filter Data
def Read_JADES(directory, Filt, Ver):
    """ Read_JADES(directory, Filt, Ver):

    Args:
        directory(string): specifies full path to filter files of specific version
        filt(string):specifies which filter we read in
        Ver(string):specifies which data version we read in

    Returns:
        image data, error data, flux density, pixel area in steradians, and pixel area in arcsec^2
    """

    #sorts files
    origDir = os.getcwd()
    os.chdir(directory)
    filelist = [f for f in os.listdir(os.getcwd()) if not f.startswith('.')]
    imagefiles={}
    if Ver== 'v1' or Ver=='v2':
        filt=str(Filt)+'_'+str(Ver)
        filename=f'hlsp_jades_jwst_nircam_goods-s-deep_{filt}.0_drz.fits'
    elif Ver=='v3':
        filt=str(Filt)+'_v1'
        filename=f'hlsp_jades_jwst_nircam_goods-n_{filt}.0_drz.fits'

    imagefiles[filt]=os.path.join(filename)

    infile = imagefiles[filt]
    print('READING JADES FILE '+ infile)
    

    global hdu
    hdu = fits.open(infile)
    
    #image data
    im = hdu[1].data
    #error ext
    err=hdu[2].data
    #define error
    err_data=err>0

    global imwcs
    imwcs = wcs.WCS(hdu[1].header, hdu)

    #primary hdu header
    hdr=hdu[0].header
    #image header
    header=hdu[1].header

    #Flux density (MJy/steradian) producing 1 cps
    photmjsr=header['PHOTMJSR']
    #Nominal pixel area in steradians
    Pix_SR=header['PIXAR_SR']
    #Nominal pixel area in arcsec^2
    Pix_A2=header['PIXAR_A2']
    #Target RA at mid time of exposure
    Targ_RA=hdr['TARG_RA']
    #Target Dec at mid time of exposure
    Targ_DEC=hdr['TARG_DEC']
    
    return im,err_data,photmjsr,Pix_SR,Pix_A2,Targ_RA,Targ_DEC


#Read in Catalog Data
def ReadJADESCatalog(dir_catalog,Ver):
    """ ReadJADESCatalog(dir_catalog, Ver):
    
    Args:
        dir_catalog(string): specifies full path to catalog files of specific version
        Ver(string):specifies which catalog version we read in

    Returns:
         ID, RA, DEC, X, Y, A(s_maj), B(s_min), PA(theta), BBOX_XMIN, BBOX_XMAX, BBOX_YMIN, BBOX_YMAX, NPIX, EXT4, EXT5, EXT6
    """
    #Read in Catalog
    if Ver=='v1' or Ver=='v2':
        hdu=fits.open(str(dir_catalog)+"\hlsp_jades_jwst_nircam_goods-s-deep_photometry_"+str(Ver)+".0_catalog.fits")
    elif Ver=='v3':
        hdu=fits.open(str(dir_catalog)+"\hlsp_jades_jwst_nircam_goods-n_photometry_v1.0_catalog.fits")
    print('READING JADES CATALOG FILE ')
    

    """EXTENSION 3     EXTNAME: SIZE
    Source parameters : ID, position, number of pixels, shape parameters, half-light radii.
    """
    #Define Ext. 3
    T3 = Table.read(hdu[3])

    #Define Table Columns for necessary information for masking 

    #• ID ‐ (slong)
    #Unique ID of source
    ID_source = T3['ID']

    #• X ‐ pix (double)
    #x‐centroid in detection image, windowed
    x=T3['X']
    #• Y ‐ pix (double)
    #y‐centroid in detection image, windowed
    y=T3['Y']

    #• A ‐ arcsec (single)
    #Semi‐major axis length, from detection image
    s_maj=T3['A']
    #• B ‐ arcsec (single)
    #Semi‐minor axis length, from detection image
    s_min=T3['B']
    #• PA ‐ deg(single)
    #position angle, from detection image
    theta= T3['PA']

    # BBOX_XMIN ‐ pix (long)  
    #Bounding box x‐minimum in detection image
    bbox_xmin=T3['BBOX_XMIN']
    #• BBOX_XMAX ‐ pix (long)
    #Bounding box x‐maximum in detection image
    bbox_xmax=T3['BBOX_XMAX']
    #• BBOX_YMIN ‐ pix (long)
    #Bounding box y‐minimum in detection image
    bbox_ymin=T3['BBOX_YMIN']
    #• BBOX_YMAX ‐ pix (long)
    #Bounding box y‐maximum in detection image
    bbox_ymax=T3['BBOX_YMAX']

    #• NPIX_DET ‐ (single)  
    #Number of pixels in the detection image belonging to the object
    npix=T3['NPIX_DET']


    """EXTENSION 4    EXTNAME: CIRC
    Circular aperture photometry, full resolution image
    """
    circ_phot= Table.read(hdu[4])
    """EXTENSION 5    EXTNAME: CIRC_BSUB
    Circular aperture photometry, background‐subtracted full resolution image
    """
    circ_bsub= Table.read(hdu[5])
    """EXTENSION 6   EXTNAME: CIRC_CONV
    Circular aperture photometry, psf‐convolved image to F444W resolution
    """
    circ_conv= Table.read(hdu[6])
    """EXTENSION 2   EXTNAME: FLAG
    """
    T2=Table.read(hdu[2])

    #• FLAG_ST ‐ Classification(integer)
    #star=1, galaxy=0
    flag_st=T2['FLAG_ST']
    
    return ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,circ_phot,circ_bsub,circ_conv,flag_st


def Read_NircamPSF(directory,Filt):
    """ Read_NircamPSF(directory,Filt):
        #Read in Nircam PSF for filter
    Args:
        directory(string): specifies full path to psf files
        Filt(string):specifies which filter we read in

    Returns:
        psf(array): psf data
    """

    #sorts files
    origDir = os.getcwd()
    os.chdir(directory)
    filelist = [f for f in os.listdir(os.getcwd()) if not f.startswith('.')]
    psf_files={}
    filt=str(Filt).upper()
    filename=f'PSF_NIRCam_in_flight_opd_filter_{filt}.fits'
    psf_files[filt]=os.path.join(filename)

    infile = psf_files[filt]
    print('READING PSF FILE '+ infile)

    hdu_psf = fits.open(infile)
    
    #psf data
    psf = hdu_psf[1].data

    return psf

def Read_TrilegalPickle(directory,Filt,Ver):
    """ Read_TrilegalPickle(directory,Filt,Ver):
        #Read in Trilegal Pickle file 
    Args:
        directory(string): specifies full path to trilegal pickle files
        Filt(string):specifies which filter we read in
        Ver(string):specifies which data version we read in

    Returns:
        trilegal_pickle(array): trilegal pickle data
    """

    #sorts files
    origDir = os.getcwd()
    os.chdir(directory)
    filelist = [f for f in os.listdir(os.getcwd()) if not f.startswith('.')]
    pickle_files={}
    filt=str(Filt).upper()
    if Ver=='v2':
        filename=f'trilegal_{filt}.pickle'
    else:
        filename=f'trilegal_{filt}_{Ver}.pickle'
    pickle_files[filt]=os.path.join(filename)

    infile = pickle_files[filt]
    print('READING TRILEGAL PICKLE FILE '+ infile)

    with open(str(filename),'rb') as file:
        trilegal_pickle=pickle.load(file)

    return trilegal_pickle

def run_Driver():
    print('---------------------------------------------------------------------------------------------------------------------------')
    
    #Read in image data
    im_data,err_data,photmjsr,Pix_SR,Pix_A2,Targ_RA,Targ_DEC=Read_JADES(Im_Dir,Filt,Ver)
    print('---------------------------------------------------------------------------------------------------------------------------')
    #Read catalog
    ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,circ_phot,circ_bsub,circ_conv,flag_st=ReadJADESCatalog(Cat_Dir,Ver)
    #Read PSF
    psf=Read_NircamPSF(PSF_Dir,Filt)
    #Read Trilegal Pickle File
    trilegal_sims=Read_TrilegalPickle(Trilegal_Pic_Dir,Filt,Ver)
    #Read vega-sirius_Jy
    zeropoint, vega_zp, sirius_vega_zp, mjsr_mean, pix_sr_mean= calc_zeropoint(Filt.upper(),file_path_zeropoints)
    print('---------------------------------------------------------------------------------------------------------------------------')
    
    #full list of filters and wavelengths for JADES
    filt_list_full=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f430m','f444w','f460m','f480m']
    wave_list_full=[.901,1.154,1.501,1.845,1.990,2.093,2.786,3.365,3.563,4.092,4.280,4.421,4.624,4.834]

    in_list=filt_list_full.index(Filt)

    #Pixel to Arcsecond Ratio
    Pix_Scale= np.sqrt(Pix_A2)
    print('Pix_Scale=',Pix_Scale)

    #define Pivot Wavelength in µm 
    Wave=wave_list_full[in_list]
    
    #Define which extension to use flux column
    #Circ aperture (0)
    if Flux_Ext==0:
        flux=circ_phot[str(Filt.upper())+'_'+str(Flux_Col)]
        flux_err=circ_phot[str(Filt.upper())+'_'+str(Flux_Col)+'_e']
    #Background-Subtracted (1)
    elif Flux_Ext==1:
        flux=circ_bsub[str(Filt.upper())+'_'+str(Flux_Col)]
        flux_err=circ_bsub[str(Filt.upper())+'_'+str(Flux_Col)+'_e']
    #Psf-convolved to f444w (2)
    elif Flux_Ext==2:
        flux=circ_conv[str(Filt.upper())+'_'+str(Flux_Col)]
        flux_err=circ_conv[str(Filt.upper())+'_'+str(Flux_Col)+'_e']

    #Define Maximum Magnitude for Scaling Eq in Masking Functions
    Max_Mag=np.ma.masked_invalid(((-2.5*np.log10(flux * (10**(-9)) ))+8.9)).max()
    print('Max_Mag=',Max_Mag)


    #Define Frequency (c/λ)
    nu= ((3*10**8)/(Wave*10**(-6)))

    #Multiply data by nosource_mask
    if Ver=='v2':
        #For these filters there are areas masked out in nosource_mask that have data but no sources in catalog
        if Filt==('f090w' or 'f115w' or 'f150w'):
            with open(r'C:\Users\kasia\nh_pipeline\nosource_mask.pickle','rb') as file:
                nosource_mask=pickle.load(file)
        else:
            nosource_mask=np.ones(np.shape(im_data))
    else:
        nosource_mask=np.ones(np.shape(im_data))
    
    print('---------------------------------------------------------------------------------------------------------------------------')

    #Flag = 1: Measures EBL for one filter
    if Flag==1:
        #set trilegal mag mask limit to 27 ab mag
        max_mask=27

        #Multiply no source mask with data
        im_data_all=nosource_mask*im_data

        #Place nan where there is no data or error
        filt_good_all=np.where(err_data,im_data_all,np.nan)

        #Calculate number of pixels
        data_area=im_data*nosource_mask
        num_pixels=np.nansum(np.isfinite(np.where(err_data,data_area,np.nan)))

        filt_good_mask=np.where(np.isnan(filt_good_all)==True,0,1)

        surface_brightness_bg= np.nanmean(((filt_good_all* 10**6)* nu * 10**(-26) * 10**9))
        print("Mean surface brightness before masking = ", surface_brightness_bg, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
    
        #Mask Image Data 
        mask,percent_pix,hist,mean_nonmask,mean_mask,mean_mask_gal_plus_dimstars=Final_Masking_new(im_data,err_data,Max_Mag,flux,flux_err,Pix_Scale,ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,0,Wave,nosource_mask,flag_st,0,max_mask)
        mask_err,err_percent_pix,hist,mean_nonmask_err,mean_mask_err,mean_mask_gal_plus_dimstars_err=Final_Masking_new(im_data,err_data,Max_Mag,flux,flux_err,Pix_Scale,ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,1,Wave,nosource_mask,flag_st,0,max_mask)
        print('---------------------------------------------------------------------------------------------------------------------------')
        
        print("MEAN mask- all")
        print("Mean surface brightness of resolved stars and galaxies = ", mean_mask, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Mean surface brightness after masking = ", mean_nonmask, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        err_resolved=np.abs(mean_mask_err-mean_mask)
        err_nonmask= np.abs(mean_nonmask_err-mean_nonmask)
        print("Error of surface brightness of resolved stars and galaxies = ", err_resolved, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Error of surface brightness after masking = ", err_nonmask, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        
        print("MEAN mask- gal")
        print("Mean surface brightness of galaxies = ", mean_mask_gal_plus_dimstars, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        # print("Mean surface brightness after masking = ", mean_nonmask_gal_plus_dimstars, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        err_resolved_gal_plus_dimstars=np.abs(mean_mask_gal_plus_dimstars_err-mean_mask_gal_plus_dimstars)
        # err_nonmask_gal_plus_dimstars= np.abs(mean_nonmask_gal_plus_dimstars_err-mean_nonmask_gal_plus_dimstars)
        print("Error of surface brightness of resolved galaxies = ", err_resolved_gal_plus_dimstars, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        # print("Error of surface brightness after masking = ", err_nonmask_gal_plus_dimstars, " nW m^-2 sr^-1", "(", Filt, Ver, ")")

        #Run Trilegal ISL to calculate surface brightness of stars greater than max_mask
        mean_tri,err_tri=jades_trilegalisl(trilegal_sims,num_pixels,Pix_A2,zeropoint,nu,max_mask)
        print("Mean surface brightness of trilegal simulations over masking limit= ", mean_tri, "nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Error of surface brightness of trilegal simulations over masking limit= ", err_tri, "nW m^-2 sr^-1", "(", Filt, Ver, ")")

        #Run input_mapper - psf wings
        sourceMap=input_mapper(len(ID_source), im_data.shape,y,x,flux,flux_err,0)
        sourceMap_err=input_mapper(len(ID_source), im_data.shape,y,x,flux,flux_err,1)

        #ISL Calculation-psf wings
        sourceMap_convolved, ISL_mean= ISL(psf,sourceMap,Wave,Pix_SR,mask)
        sourceMap_convolved, ISL_mean_err= ISL(psf,sourceMap_err,Wave,Pix_SR,mask)
        #calc psf wing error
        err_psf=np.abs(ISL_mean_err-ISL_mean)

        print("Mean surface brightness of PSF Wings (mask stars and gal)= ", ISL_mean, "nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Error of surface brightness of PSF Wings (mask stars and gal)= ", err_psf, "nW m^-2 sr^-1", "(", Filt, Ver, ")")
    
        #Mean CIB Brightness
        lam_CIB = mean_nonmask + mean_mask_gal_plus_dimstars - mean_tri - ISL_mean
        
        #Print Photometric Calibration Error
        if Filt!= 'f430m' and Filt!='f460m':
            err_photcal=.02* lam_CIB
            print("Photometric Calibration Error = ", err_photcal, "nW m^-2 sr^-1", "(", Filt, Ver, ")")
        elif Filt=='f430m' or Filt=='f460m':
            err_photcal=.04* lam_CIB
            print("Photometric Calibration Error = ", err_photcal, "nW m^-2 sr^-1", "(", Filt, Ver, ")")

        #Quadrature Error
        Quad_Error=np.sqrt((err_nonmask+ err_resolved_gal_plus_dimstars)**2 + (err_psf)**2 + (err_photcal)**2+(err_tri)**2)
        
        print("Mean CIB Brightness = ", lam_CIB, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Total CIB Error (quadrature)= ", Quad_Error, " nW m^-2 sr^-1", "(", Filt, Ver, ")")

    #Flag = 2: Measures EBL for One Chunk for a filter
    if Flag==2:
        with open(r'd:\JADES\ISL\Trilegal_pickles\v2_chunk_mask_m_11.pickle','rb') as f:
            chunk_mask=pickle.load(f)
        chunk_num=11  

        #set trilegal mag mask limit to 27 ab mag
        max_mask=27

        #Multiply no source mask with data
        im_data_all=nosource_mask*im_data

        #Place nan where there is no data or error
        filt_good_all=np.where(err_data,im_data_all,np.nan)

        #Calculate number of pixels
        data_area=im_data*nosource_mask
        num_pixels=np.nansum(np.isfinite(np.where(err_data,data_area,np.nan)))

        filt_good_mask=np.where(np.isnan(filt_good_all)==True,0,1)

        surface_brightness_bg= np.nanmean(((filt_good_all* 10**6)* nu * 10**(-26) * 10**9))
        print("Mean surface brightness before masking = ", surface_brightness_bg, " nW m^-2 sr^-1", "(", Filt, Ver, ")")

        orig_im_data=im_data.copy()
        #multiply chunk mask with data
        im_data=orig_im_data*chunk_mask
        #calculate number of pixels in chunk
        num_pixels_chunk=np.nansum(np.isfinite(np.where(err_data,im_data,np.nan)))
        #Mask Image Data 
        mask,percent_pix,hist,mean_nonmask,mean_mask,mean_mask_gal_plus_dimstars=Final_Masking_new(im_data,err_data,Max_Mag,flux,flux_err,Pix_Scale,ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,0,Wave,nosource_mask,flag_st,0,max_mask)
        mask_err,err_percent_pix,hist,mean_nonmask_err,mean_mask_err,mean_mask_gal_plus_dimstars_err=Final_Masking_new(im_data,err_data,Max_Mag,flux,flux_err,Pix_Scale,ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,1,Wave,nosource_mask,flag_st,0,max_mask)
        print('---------------------------------------------------------------------------------------------------------------------------')
        
        print("MEAN mask- all")
        print("Mean surface brightness of resolved stars and galaxies = ", mean_mask, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Mean surface brightness after masking = ", mean_nonmask, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        err_resolved=np.abs(mean_mask_err-mean_mask)
        err_nonmask= np.abs(mean_nonmask_err-mean_nonmask)
        print("Error of surface brightness of resolved stars and galaxies = ", err_resolved, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Error of surface brightness after masking = ", err_nonmask, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        
        print("MEAN mask- gal")
        print("Mean surface brightness of galaxies = ", mean_mask_gal_plus_dimstars, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        # print("Mean surface brightness after masking = ", mean_nonmask_gal_plus_dimstars, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        err_resolved_gal_plus_dimstars=np.abs(mean_mask_gal_plus_dimstars_err-mean_mask_gal_plus_dimstars)
        # err_nonmask_gal_plus_dimstars= np.abs(mean_nonmask_gal_plus_dimstars_err-mean_nonmask_gal_plus_dimstars)
        print("Error of surface brightness of resolved galaxies = ", err_resolved_gal_plus_dimstars, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        # print("Error of surface brightness after masking = ", err_nonmask_gal_plus_dimstars, " nW m^-2 sr^-1", "(", Filt, Ver, ")")

        #Run Trilegal ISL to calculate surface brightness of stars greater than max_mask
        mean_tri,err_tri=jades_trilegalisl_chunk(trilegal_sims,num_pixels,num_pixels_chunk,Pix_A2,zeropoint,nu,max_mask,chunk_num)
        print("Mean surface brightness of trilegal simulations over masking limit= ", mean_tri, "nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Error of surface brightness of trilegal simulations over masking limit= ", err_tri, "nW m^-2 sr^-1", "(", Filt, Ver, ")")

        #Run input_mapper - psf wings
        sourceMap=input_mapper(len(ID_source), im_data.shape,y,x,flux,flux_err,0)
        sourceMap_err=input_mapper(len(ID_source), im_data.shape,y,x,flux,flux_err,1)

        #ISL Calculation-psf wings
        sourceMap_convolved, ISL_mean= ISL(psf,sourceMap,Wave,Pix_SR,mask)
        sourceMap_convolved, ISL_mean_err= ISL(psf,sourceMap_err,Wave,Pix_SR,mask)
        #calc psf wing error
        err_psf=np.abs(ISL_mean_err-ISL_mean)

        print("Mean surface brightness of PSF Wings (mask stars and gal)= ", ISL_mean, "nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Error of surface brightness of PSF Wings (mask stars and gal)= ", err_psf, "nW m^-2 sr^-1", "(", Filt, Ver, ")")
    
        #Mean CIB Brightness
        lam_CIB = mean_nonmask + mean_mask_gal_plus_dimstars - mean_tri - ISL_mean
        
        #Print Photometric Calibration Error
        if Filt!= 'f430m' and Filt!='f460m':
            err_photcal=.02* lam_CIB
            print("Photometric Calibration Error = ", err_photcal, "nW m^-2 sr^-1", "(", Filt, Ver, ")")
        elif Filt=='f430m' or Filt=='f460m':
            err_photcal=.04* lam_CIB
            print("Photometric Calibration Error = ", err_photcal, "nW m^-2 sr^-1", "(", Filt, Ver, ")")

        #Quadrature Error
        Quad_Error=np.sqrt((err_nonmask+ err_resolved_gal_plus_dimstars)**2 + (err_psf)**2 + (err_photcal)**2+(err_tri)**2)
        
        print("Mean CIB Brightness = ", lam_CIB, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Total CIB Error (quadrature)= ", Quad_Error, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
    

#run_Driver()
print('------')

