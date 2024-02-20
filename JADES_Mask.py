
""" JADES_Mask.py

    Function to Mask sources in image data
    Returns mask array of 0 for sources and err and clip masking, and 1 elsewhere as first output and a list of arrays for various stages of the masking (for histogram) as second output

    1/24/24
"""

#Imports

#math
import numpy as np
import astropy.units as u
#regions-masking
from regions.core import PixCoord
from regions import EllipsePixelRegion
import numpy.ma as ma


def Final_Masking(im_data,data_err,max_mag,filt,pix_scale,ID,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix):
    
    """ Final_Masking(im_data,data_err,max_mag,filt,pix_scale,ID,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix):
        #Reads image data and error data from same filter, masks error, masks sources, runs iterations of clip masking, and returns final mask of 0s and 1s
    Args:
        im_data('numpy.ndarray'): image data
        data_err('numpy.ndarray'): error data from im ext
        max_mag('int'): maximum magnitude for masking sources (used 40 for f090w and 38 for f150w and f200w)
        filt('string'): which flux to use for AB Mag (ex. cat[13][Filt_Mag]- cat[13]= Ext.4, cat[14]= Ext.5 (bsub), cat[15]=Ext.6 (conv)
                                                        [Filt_Mag] specified in JADES_config ex. 'F090W_CIRC0')
        pix_scale('float'): Arcsecond to Pixel Ratio (0.031"/pix for 0.6-2.3 µm and 0.063"/pix for 2.4-5 µm )
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
    Returns:
        fin_mask (class 'numpy.ndarray'): Returns mask array with 0 for sources and err and clip masking, and 1 elsewhere
        listhist ('list'): list of arrays to plot in histogram including before masking (excluding error), after masking sources, and after multiple rounds of clip masks
    """  
    
    #create list to append arrays of different stages of masking for histogram
    listhist=[]
    #Place nan where there is no data or error
    filt_good=np.where(data_err,im_data,np.nan)
    #Append array of image data before any masking to show in histogram
    listhist.append(filt_good)

    #create mask where error becomes 0 and otherwise is 1
    obj_mask=np.where(data_err,1,0)

    #calculate AB magnitude of flux
    flux=(-2.5*np.log10(filt * (10**(-9)) ))+8.9
    
    #over the range of the catalog (all the sources)
    for i in range(len(ID)):
        
        #retrieve the magnitude (Jy) of flux from specified filter for each source
        thismag=flux[i]

        #if the magnitude is a valid number (not infinite or nan) continue with filter scaling
        if thismag != np.inf and np.isnan(thismag)!=True:
            #calculate scaling factor
            scaling=2.5* (max_mag/ thismag)**2
        
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
            #determine for each source whether each pixel in im_source shape are inside or outisde ellipse pixel region   
            object_mask_contains = object_aperture.contains(PixCoord(x=np.arange(im_source.shape[1]),y=np.arange(im_source.shape[0])[:,np.newaxis]))
            #if pixel is inside region it is assigned as 0
            obj_mask_source[object_mask_contains]=0
            #puts obj mask array of source into main obj mask array for image by mulitiplying each pixel by 0 if inside region and 1 if outside in order to preserve zeros
            obj_mask[ymin_source:ymax_source, xmin_source:xmax_source]=obj_mask_source*obj_mask[ymin_source:ymax_source, xmin_source:xmax_source]

    
    #create new array = to im.data everywhere except where obj_mask=0, in which case those values are assigned np.nan in order to avoid forced zeros in clip masking
    im_new=np.where(obj_mask==0,np.nan,im_data)
    #redefine im_new as masked array in order to mask invalid values (nan and infs) to exclude these values from calculating mean or standard deviation
    im_new=ma.masked_invalid(im_new)
    #Append data array after source masking to show in histogram
    listhist.append(im_new.copy())
    
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
    
    #return final mask where all masked values become 0 and elsewhere is 1
    fin_mask=np.where(im_new.mask==True,0,1)

    #print percentage of original pixels left
    print("The percent of original pixels left = ",(np.sum(np.isfinite(im_new))/np.sum(np.isfinite(filt_good)))*100)
    return fin_mask, listhist