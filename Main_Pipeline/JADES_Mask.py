
""" JADES_Mask.py

    Function to Mask sources and calculate surface brightness

    9/28/24
"""

#Imports

#math
import numpy as np
import astropy.units as u
import pickle
#regions-masking
from regions.core import PixCoord
from regions import EllipsePixelRegion
import numpy.ma as ma


def Final_Masking_new(im_data,data_err,max_mag,flux_cat,fluxerr,pix_scale,ID,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,err_flag,wavelength,nosource_mask, star_flag, no_negpix,max_mask):
    """ Final_Masking_new():
        #Reads image data and error data from a filter, masks error, masks sources, runs iterations of clip masking, returns final mask of 0s and 1s, and calculates surface brightness of components
    Args:
        im_data('numpy.ndarray'): image data
        data_err('numpy.ndarray'): error data
        max_mag('int'): maximum magnitude for masking sources
        flux_cat('string'): flux column from catalog
        fluxerr('string'): flux uncertainty column
        pix_scale('float'): Arcsecond to Pixel Ratio
        ID('numpy.ndarray'): • ID - (slong), Unique ID of source
        x('numpy.ndarray'): • X - pix (double), x-centroid in detection image, windowed
        y('numpy.ndarray'): • Y - pix (double), y-centroid in detection image, windowed
        s_maj('numpy.ndarray'): • A - arcsec (single), Semi-major axis length, from detection image
        s_min('numpy.ndarray'): • B - arcsec (single), Semi-minor axis length, from detection image
        theta('numpy.ndarray'): • PA - deg(single), position angle, from detection image
        bbox_xmin('numpy.ndarray'): • BBOX_XMIN - pix (long), Bounding box x-minimum in detection image
        bbox_xmax('numpy.ndarray'): • BBOX_XMAX - pix (long), Bounding box x-maximum in detection image
        bbox_ymin('numpy.ndarray'): • BBOX_YMIN - pix (long), Bounding box y-minimum in detection image
        bbox_ymax('numpy.ndarray'): • BBOX_YMAX - pix (long), Bounding box y-maximum in detection image
        npix('numpy.ndarray'): • NPIX_DET - (single), Number of pixels in the detection image belonging to the object
        err_flag(0 or 1): 0= regular flux calculation, 1= include error calculation
        wavelength('float'): Pivot Wavelength in micrometers
        nosource_mask('numpy.ndarray'): mask of areas with data but no sources
        star_flag('numpy.ndarray'): classification of stars or galaxies in catalog)
        no_negpix(0 or 1): 1- get rid of negative pixels for mean, 0-keep negative pixels
        max_mask('int'): include stars with a magnitude greater than this number
    Returns:
        fin_mask('numpy.ndarray'): Returns mask array with 0 for stars and galaxies and err and clip masking, and 1 elsewhere
        percent_pix('float'): percent of original pixels left after all sources are masked + clip masking
        listhist('list'): list of arrays to plot in histogram including before masking (excluding error), after masking sources, and after multiple rounds of clip masks
        surface_brightness_bg('float'): Surface brightness of background when all sources are masked + clip mask
        surface_brightness_mask('float'): Surface brightness inside mask of all sources
        surface_brightness_gal_plus_dimstars('float'): Surface brightness inside mask of galaxies and dim stars with magnitude greater than max_mask
    """  
    #seed random numbers used for error calculation
    np.random.seed(1)

    #create list to append arrays of different stages of masking for histogram
    listhist=[]
    
    #Multiply no source mask with data
    im_data=nosource_mask*im_data

    #Place nan where there is no data or error
    filt_good=np.where(data_err,im_data,np.nan)
    if no_negpix==1:
        filt_good=np.where(filt_good<0,np.nan,filt_good)
        
    im_data= filt_good.copy()

    #Append array of image data before any masking to show in histogram
    listhist.append(filt_good)

    #create masks where error becomes 0 and otherwise is 1
    obj_mask=np.where(data_err,1,0)
    obj_mask_gal_plus_dimstars=obj_mask.copy()
    
    #over the range of the catalog (all the sources)
    for i in range(len(ID)):
        
        if err_flag==0:
            #retrieve the magnitude (nJy) of flux from specified filter for each source
            flux=flux_cat[i]
        if err_flag==1:

            thisfluxerr=fluxerr[i]
            #if masked/no error value
            if ma.is_masked(thisfluxerr)==True:
                thisfluxerr=0
            #retrieve the magnitude (nJy) of flux from specified filter for each source
            flux=flux_cat[i]+ (thisfluxerr * np.random.randn(1)[0])
        
        #calculate AB magnitude of flux
        thismag=(-2.5*np.log10(flux * (10**(-9)) ))+8.9

        #if the magnitude is a valid number (not infinite or nan) continue with filter scaling
        if np.abs(thismag) != np.inf and np.isnan(thismag)!=True:
            
            #calculate scaling factor
            scaling=2.5* (max_mag/ thismag)**2

        #if magnitude is inf or nan continue with npix scaling using catalog     
        else:
            
            #calculate area of each source as provided in catalog with npix (number of pixels in the detection image belonging to object)
            #area=(s_min/pix_scale)*scaling*(s_maj/pix_scale)*scaling
            area=npix[i]
            #calculate scale of major and minor axis
            s_maj_scale=s_maj[i]/pix_scale
            s_min_scale=s_min[i]/pix_scale
            mult=s_maj_scale*s_min_scale
            #calculate scaling factor
            scaling=(area/mult)**(1/2)
            
        #define a box for each object [i] based on bounding box +- 50 in order to shorten area that object_mask_contains searches if each pixel is inside or outside region
        xmin_source = bbox_xmin[i]-50
        xmax_source = bbox_xmax[i]+50
        ymin_source = bbox_ymin[i]-50
        ymax_source = bbox_ymax[i]+50

        #cut the im array to the limits of each object
        im_source = im_data[ymin_source:ymax_source, xmin_source:xmax_source]
        #also create obj mask array for source of the same shape as im array cut to source filled with 1
        obj_mask_source=np.ones(im_source.shape)
        
        #define region of each source
        #define center of source in pixel coords
        object_center_pix = PixCoord(x[i]-xmin_source, y[i]-ymin_source)
        #create ellipse pixel region using center, height and width using major and minor axis multiplied with scaling factor calculated (and divided by .031 to translate to pixels from arcsecond), and additionally the angle of ellipse
        object_aperture=EllipsePixelRegion(object_center_pix,height=scaling*s_min[i]/pix_scale,width=scaling*s_maj[i]/pix_scale,angle=theta[i]*u.deg,visual={'color': 'white'}) 
        #determine for each source whether each pixel in im_source is inside or outisde ellipse pixel region   
        object_mask_contains = object_aperture.contains(PixCoord(x=np.arange(im_source.shape[1]),y=np.arange(im_source.shape[0])[:,np.newaxis]))
        #if pixel is inside ellipse region it is assigned as 0
        obj_mask_source[object_mask_contains]=0
        #puts obj mask array of source into main obj mask array for image by mulitiplying each pixel by 0 if inside region and 1 if outside in order to preserve zeros
        obj_mask[ymin_source:ymax_source, xmin_source:xmax_source]=obj_mask_source*obj_mask[ymin_source:ymax_source, xmin_source:xmax_source]
        
        #put obj mask array of source into main obj mask array for galaxies and dim stars
        if star_flag[i]==1:
            #star
            #If star of valid magnitude is greater than max_mask, add to gal mask 
            if np.isfinite(thismag)==True:
                if thismag>max_mask:
                    obj_mask_gal_plus_dimstars[ymin_source:ymax_source, xmin_source:xmax_source]=obj_mask_source*obj_mask_gal_plus_dimstars[ymin_source:ymax_source, xmin_source:xmax_source]
        else:
            #gal
            obj_mask_gal_plus_dimstars[ymin_source:ymax_source, xmin_source:xmax_source]=obj_mask_source*obj_mask_gal_plus_dimstars[ymin_source:ymax_source, xmin_source:xmax_source]

    #create new array = to im.data everywhere except where obj_mask=0, in which case those values are assigned np.nan in order to avoid forced zeros in clip masking
    im_new=np.where(obj_mask==0,np.nan,im_data)
    #redefine im_new as masked array in order to mask invalid values (nan and infs) to exclude these values from calculating mean or standard deviation
    im_new=ma.masked_invalid(im_new)
    #Append data array after source masking to show in histogram
    listhist.append(im_new.copy())

    #For gal- create new array = to im.data everywhere except where obj_mask=0, in which case those values are assigned np.nan
    im_new_gal_plus_dimstars=np.where(obj_mask_gal_plus_dimstars==0,np.nan,im_data)
    #For gal- redefine im_new as masked array in order to mask invalid values (nan and infs) to exclude these values from calculating mean or standard deviation
    im_new_gal_plus_dimstars=ma.masked_invalid(im_new_gal_plus_dimstars)

    #Begin clip masking
    #set counter at 0 to start
    counter=0
    #set the number of standard deviations as 3
    num_std=3

    #Run clip masking for 2 iterations
    while counter<=2:
        #calculate mean of im_new, use np.nanmean to exclude masked values
        im_mean=np.nanmean(im_new)
        #calculate standard deviation of im_new, use np.nanstd to exclude masked values
        im_std=np.nanstd(im_new)
        #where the absolute value of im_new is greater than num_std times the standard deviation plus the mean then replace value with nan
        im_new=np.where(np.absolute(im_new)<(im_mean+(num_std*im_std)),im_new,np.nan)
        #mask invalid numbers (nan), redefine im_new for next iteration
        im_new=ma.masked_invalid(im_new)
        
        #add 1 to counter to count iteration
        counter += 1
        #Append data array after each clip mask to show in histogram
        listhist.append(im_new.copy())
    
    #surface brightness of image after masking
    surface_brightness_bg= np.nanmean(((im_new* 10**6)* ((3*10**8)/(wavelength*10**(-6)))* 10**(-26) * 10**9))

    #return final mask where all masked values become 0 and elsewhere is 1
    fin_mask=np.where(im_new.mask==True,0,1)

    #opposite mask of np.nan of non-masked values, and elsewhere is original values after masking error data
    opp_mask=np.where(im_new.mask==False,np.nan,filt_good)
    #surface brightness of resolved stars and galaxies masked
    surface_brightness_mask= np.nanmean(((opp_mask* 10**6)* ((3*10**8)/(wavelength*10**(-6)))* 10**(-26) * 10**9))

    #print percentage of original pixels left
    percent_pix=(np.sum(np.isfinite(im_new))/np.sum(np.isfinite(filt_good)))*100
    print("The percent of original pixels left = ",percent_pix)

    #opposite mask of np.nan of non-masked values, and elsewhere is original values after masking error data
    opp_mask_gal_plus_dimstars=np.where(im_new_gal_plus_dimstars.mask==False,np.nan,filt_good)
    #surface brightness of galaxies
    surface_brightness_gal_plus_dimstars= np.nanmean(((opp_mask_gal_plus_dimstars* 10**6)* ((3*10**8)/(wavelength*10**(-6)))* 10**(-26) * 10**9)) 

    return fin_mask, percent_pix, listhist, surface_brightness_bg, surface_brightness_mask,surface_brightness_gal_plus_dimstars