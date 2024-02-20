
""" JADES_driver.py

    Reads in JADES image data for filter, and catalog file for a version of data
    Use JADES_config file in order to set variables to use in this file
        Set Flag= (0- nothing, 1- Zscale with LinearStretch plot of filter image data, 2- Plot of Mask, 3- Plot of Mask*Data, 4- Histogram, 5-Chunking Plot of Data, 6- Chunking Plot of Mask*Data, 7 Chunk Histogram, 8 ISL)

    2/19/24
"""

#Imports
from astropy.io import fits
import os
from astropy.table import Table
import astropy.wcs as wcs
from astropy.io import fits
import numpy as np
#plotting
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval,LinearStretch,ImageNormalize


#Read in filter Data
def Read_JADES(directory, Filt):
    """ Read_JADES(directory, filt):

    Args:
        directory(string): specifies full path to filter files of specific version
        filt(string):specifies which filter we read in + data version

    Returns:
        output_Read_JADES (list): list containing image data, error data, imwcs data, flux density, pixel area in steradians, and pixel area in arcsec^2
    """

    #sorts files
    origDir = os.getcwd()
    os.chdir(directory)
    filelist = [f for f in os.listdir(os.getcwd()) if not f.startswith('.')]
    imagefiles={}
    filt=str(Filt)+'_'+str(Ver)
    filename=f'hlsp_jades_jwst_nircam_goods-s-deep_{filt}.0_drz.fits'
    imagefiles[filt]=os.path.join(filename)

    infile = imagefiles[filt]
    print('READING JADES FILE '+ infile)
    
    #List for outputs
    output_Read_JADES=[]

    global hdu
    hdu = fits.open(infile)
    
    #image data
    im = hdu[1]
    output_Read_JADES.append(im)
    
    #error ext
    err=hdu[2]
    output_Read_JADES.append(err)

    global imwcs
    imwcs = wcs.WCS(hdu[1].header, hdu)
    output_Read_JADES.append(imwcs)

    #image header
    header=im.header

    #Flux density (MJy/steradian) producing 1 cps
    PHOTMJSR=header['PHOTMJSR']
    output_Read_JADES.append(PHOTMJSR)
    #aNominal pixel area in steradians
    Pix_SR=header['PIXAR_SR']
    output_Read_JADES.append(Pix_SR)
    #Nominal pixel area in arcsec^2
    Pix_A2=header['PIXAR_A2']
    output_Read_JADES.append(Pix_A2)
    
    return output_Read_JADES

#Read in Catalog Data
def ReadJADESCatalog(dir_catalog, ver):
    """ ReadJADESCatalog(dir_catalog, ver):

    Args:
        dir_catalog(string): specifies full path to catalog files of specific version
        ver(string):specifies data version ('v1' or 'v2')

    Returns:
        output_ReadJADESCatalog (list): list containing ID, RA, DEC, X, Y, A(S_maj), B(S_min), 
                                        PA, BBOX_XMIN, BBOX_XMAX, BBOX_YMIN, BBOX_YMAX, NPIX, 
                                        EXT4, EXT5, EXT6
    """
    #Read in Catalog
    hdu=fits.open(str(dir_catalog)+"\hlsp_jades_jwst_nircam_goods-s-deep_photometry_"+str(ver)+".0_catalog.fits")
    print('READING JADES CATALOG FILE ')
    
    #List for outputs for Catalog
    output_ReadJADESCatalog=[]

    """EXTENSION 3     EXTNAME: SIZE
    Source parameters : ID, position, number of pixels, shape parameters, half-light radii.
    """
    #Define Ext. 3
    T3 = Table.read(hdu[3])

    #Define Table Columns for necessary information for masking 

    #• ID ‐ (slong)
    #Unique ID of source
    ID = T3['ID']
    output_ReadJADESCatalog.append(ID)

    #• RA ‐ degree (double)
    #RA of source, ICRS, from detection image
    ra = T3['RA']
    output_ReadJADESCatalog.append(ra)

    #• DEC ‐ degree (double)
    #Dec of source, ICRS, from detection image
    dec = T3['DEC']
    output_ReadJADESCatalog.append(dec)

    #• X ‐ pix (double)
    #x‐centroid in detection image, windowed
    x=T3['X']
    output_ReadJADESCatalog.append(x)

    #• Y ‐ pix (double)
    #y‐centroid in detection image, windowed
    y=T3['Y']
    output_ReadJADESCatalog.append(y)

    #• A ‐ arcsec (single)
    #Semi‐major axis length, from detection image
    s_maj=T3['A']
    output_ReadJADESCatalog.append(s_maj)

    #• B ‐ arcsec (single)
    #Semi‐minor axis length, from detection image
    s_min=T3['B']
    output_ReadJADESCatalog.append(s_min)

    #• PA ‐ deg(single)
    #position angle, from detection image
    theta= T3['PA']
    output_ReadJADESCatalog.append(theta)

    # BBOX_XMIN ‐ pix (long)  
    #Bounding box x‐minimum in detection image
    bbox_xmin=T3['BBOX_XMIN']
    output_ReadJADESCatalog.append(bbox_xmin)

    #• BBOX_XMAX ‐ pix (long)
    #Bounding box x‐maximum in detection image
    bbox_xmax=T3['BBOX_XMAX']
    output_ReadJADESCatalog.append(bbox_xmax)

    #• BBOX_YMIN ‐ pix (long)
    #Bounding box y‐minimum in detection image
    bbox_ymin=T3['BBOX_YMIN']
    output_ReadJADESCatalog.append(bbox_ymin)

    #• BBOX_YMAX ‐ pix (long)
    #Bounding box y‐maximum in detection image
    bbox_ymax=T3['BBOX_YMAX']
    output_ReadJADESCatalog.append(bbox_ymax)

    #• NPIX_DET ‐ (single)  
    #Number of pixels in the detection image belonging to the object
    npix=T3['NPIX_DET']
    output_ReadJADESCatalog.append(npix)


    """EXTENSION 4    EXTNAME: CIRC
    Circular aperture photometry, full resolution image
    """
    T4= Table.read(hdu[4])
    output_ReadJADESCatalog.append(T4)

    """EXTENSION 5    EXTNAME: CIRC_BSUB
    Circular aperture photometry, background‐subtracted full resolution image
    """
    T5= Table.read(hdu[5])
    output_ReadJADESCatalog.append(T5)

    """EXTENSION 6   EXTNAME: CIRC_CONV
    Circular aperture photometry, psf‐convolved image to F444W resolution
    """
    T6= Table.read(hdu[6])
    output_ReadJADESCatalog.append(T6)
    
    return output_ReadJADESCatalog

#Import configuration variables used to read in image data, catalog data, masking, and histogram
from JADES_config import *

#Read in image data and catalog data
dat=Read_JADES(Im_Dir,Filt)
cat=ReadJADESCatalog(Cat_Dir, Ver)


#list of filters to assign pixel to arcsecond ratio
filt_list=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f430m','f444w','f460m','f480m']
#Find index of filter given
in_list=filt_list.index(Filt)
#Pixel to Arcsecond Ratio (0.031"/pix for 0.6-2.3 µm and 0.063"/pix for 2.4-5 µm )
""" 
    0.031"/pix for F090W, F115W, F150W, F182M, F200W, F210M
    0.063"/pix for F277W, F335M, F356W, F410M, F430M, F444W, F460M, F480M
"""
#Assign pixel to arcsecond ratio for correct filter
if in_list<= 5:
    Pix_Scale=0.031
else:
    Pix_Scale=0.063

#Pivot Wavelength in µm
""" 
    F090W: 0.901, F115W: 1.154, F150W: 1.501, F182M: 1.845, F200W: 1.990, F210M: 2.093
    F277W: 2.786, F335M: 3.365, F356W: 3.563, F410M: 4.092, F430M: 4.280, F444W: 4.421, F460M: 4.624, F480M: 4.834
"""
#list of wavelength values to assign wavelength to index of filters
wave_list=[.901,1.154,1.501,1.845,1.990,2.093,2.786,3.365,3.563,4.092,4.280,4.421,4.624,4.834]
#Assign wavelength value corresponding to index of filter in list
Wave=wave_list[in_list]


#For Flag=1 Plot Image data in Zscale with LinearStretch 
if Flag==1:
    fig=plt.figure(figsize=(12, 12))
    ax=plt.subplot(projection=imwcs)
    filt_good=np.where(dat[1].data>0,dat[0].data,np.nan)
    cmap=plt.cm.viridis
    norm = ImageNormalize(filt_good, interval=ZScaleInterval(),stretch=LinearStretch())
    plt.imshow(filt_good, norm=norm, cmap=cmap,interpolation='none',origin='lower')

    plt.ylabel('Declination',fontsize=30)
    plt.xlabel('Right Ascension',fontsize=30)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15)
    if save_fig== 'Y':
        plt.savefig(file_save,dpi=dpi_num)
    plt.show()

#import Masking Function
from JADES_Mask import *

#Define what extension for flux values (T4: Circ aperture (0), T5: Bsubb- background subtracted (1), T6: Conv- psf-convolved to f444w resolution(2), default 0)
if Flux_Col==0:
    output_flux=13
elif Flux_Col==1:
    output_flux=14
elif Flux_Col==2:
    output_flux=15

#Define Maximum Magnitude for Scaling Eq in Masking Functions
Max_Mag=np.ma.masked_invalid(((-2.5*np.log10(cat[output_flux][str(Filt.upper())+'_'+str(Flux)] * (10**(-9)) ))+8.9)).max()

#For Flag=2 Plot mask of 0s and 1s
if Flag==2:
    #Mask Image Data with inputs of im_data,data_err,max_mag,filt,pix_scale,ID,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix
    mask=Final_Masking(dat[0].data,dat[1].data>0,Max_Mag,cat[output_flux][str(Filt.upper())+'_'+str(Flux)],Pix_Scale,cat[0],cat[3],cat[4],cat[5],cat[6],cat[7],cat[8],cat[9],cat[10],cat[11],cat[12])[0]
    
    fig = plt.figure(figsize=(12, 12)) 
    ax1 = plt.subplot(projection=imwcs)
    plt.imshow(mask, cmap=plt.cm.viridis,interpolation='nearest', origin='lower')

    plt.ylabel('Declination',fontsize=30)
    plt.xlabel('Right Ascension',fontsize=30)
    plt.rcParams.update({'font.size': 20})
    fig.tight_layout()
    plt.subplots_adjust(left=0.15)

    if save_fig== 'Y':
        plt.savefig(file_save,dpi=dpi_num)
    plt.show()

#For Flag=3 Plot mask*data
if Flag==3:
    #Create mask with inputs of im_data,data_err,max_mag,filt,pixscale,ID,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix
    mask=Final_Masking(dat[0].data,dat[1].data>0,Max_Mag,cat[output_flux][str(Filt.upper())+'_'+str(Flux)],Pix_Scale,cat[0],cat[3],cat[4],cat[5],cat[6],cat[7],cat[8],cat[9],cat[10],cat[11],cat[12])[0]
    
    fig=plt.figure(figsize=(12, 12))
    ax=plt.subplot(projection=imwcs)
    #Where masked place nan
    new_arr=np.where(mask==0,np.nan,mask)
    norm = ImageNormalize(hdu[1].data*new_arr, interval=ZScaleInterval(),vmin=0,stretch=LinearStretch())
    plt.imshow(hdu[1].data*new_arr, cmap='gray', norm=norm, interpolation='none',origin='lower')

    plt.ylabel('Declination',fontsize=30)
    plt.xlabel('Right Ascension',fontsize=30)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15)

    if save_fig== 'Y':
        plt.savefig(file_save,dpi=dpi_num)
    plt.show()

#For Flag=4 Create Histogram
if Flag==4:
    #Create histogram with inputs of im_data,data_err,max_mag,filt,pixscale,ID,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix
    hist=Final_Masking(dat[0].data,dat[1].data>0,Max_Mag,cat[output_flux][str(Filt.upper())+'_'+str(Flux)],Pix_Scale,cat[0],cat[3],cat[4],cat[5],cat[6],cat[7],cat[8],cat[9],cat[10],cat[11],cat[12])[1]
    
    #Calculate frequency (c (m/s)/ wavelength (m))                                                                            
    nu=((3*10**8)/(Wave*10**(-6)))
    #Create list to translate array from histograms to correct units
    h_list=[]
    #For each array translate to units of nW/ (m^2 sr)
    for val in hist:
        h_list.append(val*nu*10**(-11))
    #Plot Histogram
    fig=plt.figure(figsize=(10,10))
    plt.hist(h_list[0].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,label='Image data')
    plt.hist(h_list[1].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,color='turquoise',label='After sources masked')
    plt.hist(h_list[2].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,color='salmon',label='3σ clip mask')
    plt.hist(h_list[3].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,color='yellow',label='Clip mask 2nd iteration')
    plt.hist(h_list[4].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,color='darkorchid',label='Clip mask 3rd iteration')

    plt.xlabel('Image Data [nW/(m^2 sr)]',fontsize=16)
    plt.ylabel('Count',fontsize=16)
    plt.title(str(Filt.upper())+'_'+str(Flux)+'_'+str(Ver))
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.rcParams.update({'font.size': 20})
    plt.legend()

    if save_fig== 'Y':
        plt.savefig(file_save,dpi=dpi_num)
    plt.show()

from JADES_Chunking import *
#For Flag=5 Chunking Plot of just Box
if Flag==5:
    fig=plt.figure(figsize=(12, 12))
    ax=plt.subplot(projection=imwcs)
    filt_good=np.where((dat[1].data>0)[y_min:y_max, x_min:x_max],dat[0].data[y_min:y_max, x_min:x_max],np.nan)
    #cmap=plt.cm.viridis
    norm = ImageNormalize(filt_good, interval=ZScaleInterval(),stretch=LinearStretch())
    plt.imshow(filt_good, norm=norm, cmap='gray',interpolation='none',origin='lower')

    plt.ylabel('Declination',fontsize=30)
    plt.xlabel('Right Ascension',fontsize=30)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15)
    if save_fig== 'Y':
        plt.savefig(file_save,dpi=dpi_num)
    plt.show()

#For Flag=6 Chunking Plot of Mask*data
if Flag==6:
    #Create mask with inputs of im_data,data_err,max_mag,filt,pixscale,ID,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix
    mask=Final_Masking_Chunk(x_min,x_max,y_min,y_max,dat[0].data,dat[1].data>0,Max_Mag,cat[output_flux][str(Filt.upper())+'_'+str(Flux)],Pix_Scale,cat[0],cat[3],cat[4],cat[5],cat[6],cat[7],cat[8],cat[9],cat[10],cat[11],cat[12])[0]
    
    fig=plt.figure(figsize=(12, 12))
    ax=plt.subplot(projection=imwcs)
    #Where masked place nan
    new_arr=np.where(mask==0,np.nan,mask)
    cmap=plt.cm.viridis
    cmap.set_bad(color='white')
    norm = ImageNormalize(dat[0].data[y_min:y_max, x_min:x_max]*new_arr, interval=ZScaleInterval(),vmin=0,stretch=LinearStretch())
    plt.imshow(dat[0].data[y_min:y_max, x_min:x_max]*new_arr, cmap=cmap, norm=norm, interpolation='none',origin='lower')

    plt.ylabel('Declination',fontsize=30)
    plt.xlabel('Right Ascension',fontsize=30)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15)

    if save_fig== 'Y':
        plt.savefig(file_save,dpi=dpi_num)
    plt.show()

#For Flag=7 Create Chunking Histogram
if Flag==7:
    #Create histogram with inputs of im_data,data_err,max_mag,filt,pixscale,ID,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix
    hist=Final_Masking_Chunk(x_min,x_max,y_min,y_max,dat[0].data,dat[1].data>0,Max_Mag,cat[output_flux][str(Filt.upper())+'_'+str(Flux)],Pix_Scale,cat[0],cat[3],cat[4],cat[5],cat[6],cat[7],cat[8],cat[9],cat[10],cat[11],cat[12])[1]
    
    #Calculate frequency (c (m/s)/ wavelength (m))                                                                            
    nu=((3*10**8)/(Wave*10**(-6)))
    #Create list to translate array from histograms to correct units
    h_list=[]
    #For each array translate to units of nW/ (m^2 sr)
    for val in hist:
        h_list.append(val*nu*10**(-11))

    #Plot Histogram
    fig=plt.figure(figsize=(10,10))
    plt.hist(h_list[0].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,label='Image data')
    plt.hist(h_list[1].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,color='turquoise',label='After sources masked')
    plt.hist(h_list[2].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,color='salmon',label='3σ clip mask')
    plt.hist(h_list[3].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,color='yellow',label='Clip mask 2nd iteration')
    plt.hist(h_list[4].flatten(),bins=1000,range=[-.05*nu*10**(-11),.1*nu*10**(-11)],log=True,color='darkorchid',label='Clip mask 3rd iteration')

    plt.xlabel('Image Data [nW/(m^2 sr)]',fontsize=16)
    plt.ylabel('Count',fontsize=16)
    plt.title(str(Filt.upper())+'_'+str(Flux)+'_'+str(Ver)+': '+'Chunk: '+'('+str(y_min)+','+str(y_max)+':'+str(x_min)+','+str(x_max)+')')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.rcParams.update({'font.size': 20})
    plt.legend()

    if save_fig== 'Y':
        plt.savefig(file_save,dpi=dpi_num)
    plt.show()
    
#For Flag=8 ISL
if Flag==8:
    print('PHOTMJSR=',dat[3], ', PIX_SR=',dat[4],', PIX_A2=',dat[5])
    
print('------')

