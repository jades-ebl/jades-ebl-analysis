""" Make_Im.py

    define functions to sum intensity per magnitude bin for sources
    run with Make_Im_Run.py

    8/26/24
"""

#Imports

#math
import numpy as np
import astropy.units as u

#regions-masking
from regions.core import PixCoord
from regions import EllipsePixelRegion


#Sum Flux of sources by masking a noiseless simulation and calculate surface brightness
def Sum_Flux_Mask(omega_pix,sourceMap_convolved,im_data,data_err,max_mag,flux_cat,fluxerr,pix_scale,ID,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,wavelength,nosource_mask,zpvega,mjsr,pix_sr,star_col,star_or_gal):    
    """ Sum_Flux_Mask():

    Args:
        omega_pix('float'): Pixel Area in Steradians
        sourceMap_convolved('numpy.ndarray'): noiseless simulation
        im_data('numpy.ndarray'): image data
        data_err('numpy.ndarray'): error data from im ext
        max_mag('int'): maximum magnitude in catalog
        flux_cat('numpy.ndarray'): flux column from catalog
        fluxerr('numpy.ndarray'): corresponded flux error from catalog
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
        wavelength('float'): Pivot Wavelength in micrometers
        nosource_mask('numpy.ndarray'): mask of areas with data but no sources
        zpvega('float'): vega zeropoint
        mjsr('float'): Flux density (MJy/steradian) producing 1 cps
        pix_sr('float'): pixel area in steradians
        star_col('numpy.ndarray'): classification of stars or galaxies in catalog
        star_or_gal(0 or 1): 0-only stars, 1-only galaxies

    Returns:
        sources_list_cat('list'): list of AB magnitudes of sources
        sources_list_cat_vega('list'): list of vega magnitudes of sources
        sb_list('list'): list of surface brightness of sources in nW m^-2 sr^-1
        err_list('list'): list of error of surface brightness of sources in nW m^-2 sr^-1
    """

    #create lists to append magnitudes and surface brightness
    sources_list_cat=[]
    sources_list_cat_vega=[]
    sb_list=[]
    err_list=[]

    #Multiply no source mask with data
    im_data=nosource_mask*im_data

    #create mask where error becomes np.nan and otherwise is 1
    obj_mask=np.where(data_err,1,np.nan)
    obj_mask=nosource_mask*obj_mask

    #calculate AB magnitude of flux, and vega magnitudes
    flux=(-2.5*np.log10(flux_cat * (10**(-9)) ))+8.9
    #vega flux in nJy
    flux_vega=((mjsr*pix_sr)/ 10**(zpvega/(-2.5)))*(10**15)
    mag_vega=-2.5*np.log10(flux_cat/flux_vega)
    
    #over the range of the catalog
    for i in range(len(ID)):
        thismag=flux[i]
        thiserr=fluxerr[i]

        #if the magnitude is a valid number (not infinite or nan) continue with filter scaling
        if thismag != np.inf and np.isnan(thismag)!=True:
            #calculate scaling factor
            scaling=2.5* (max_mag/ thismag)**2

        #if magnitude is inf or nan continue with npix scaling using catalog     
        else:
            #calculate area of each source as provided in catalog with npix (number of pixels in the detection image belonging to object)
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
        #also create obj mask array for source of the same shape as im array cut to source filled with 0s
        obj_mask_source=np.zeros(im_source.shape)

        #define region of each source
        #define center of source in pixel coords
        object_center_pix = PixCoord(x[i]-xmin_source, y[i]-ymin_source)
        #create ellipse pixel region using center, height and width using major and minor axis multiplied with scaling factor calculated (and divided by .031 to translate to pixels from arcsecond), and additionally the angle of ellipse
        object_aperture=EllipsePixelRegion(object_center_pix,height=scaling*s_min[i]/pix_scale,width=scaling*s_maj[i]/pix_scale,angle=theta[i]*u.deg,visual={'color': 'white'}) 
        #determine for each source whether each pixel in im_source is inside or outisde ellipse pixel region   
        object_mask_contains = object_aperture.contains(PixCoord(x=np.arange(im_source.shape[1]),y=np.arange(im_source.shape[0])[:,np.newaxis]))
        #if pixel is inside ellipse region it is assigned as 1
        obj_mask_source[object_mask_contains]=1
        
        #sum flux in region defined of noiseless simulation
        source_flux=np.nansum(sourceMap_convolved[ymin_source:ymax_source, xmin_source:xmax_source]*obj_mask_source)
        #calculate surface brightness of source and its error in nW m^-2 sr^-1
        sb=((source_flux* 10**-9)* ((3*10**8)/(wavelength*10**(-6)))* 10**(-26) * 10**9)/omega_pix
        err=((thiserr* 10**-9)* ((3*10**8)/(wavelength*10**(-6)))* 10**(-26) * 10**9)/omega_pix
        
        #Append values of AB magnitude, vega magnitude, surface brightness, and error
        #only stars
        if star_or_gal==0:
            
            if star_col[i]==1:
                sources_list_cat.append(flux[i])
                sources_list_cat_vega.append(mag_vega[i])
                sb_list.append(sb)
                err_list.append(err)

        #only galaxies
        if star_or_gal==1:
            
            if star_col[i]!=1:
                sources_list_cat.append(flux[i])
                sources_list_cat_vega.append(mag_vega[i])
                sb_list.append(sb)
                err_list.append(err)
    
    return sources_list_cat,sources_list_cat_vega,sb_list,err_list


#Flux of sources from catalog and calculate surface brightness
def Sum_Flux_Mask_cat(omega_pix,flux_cat,fluxerr,ID,wavelength,zpvega,mjsr,pix_sr,star_col,star_or_gal):    
    """ Sum_Flux_Mask_cat():

    Args:
        omega_pix('float'): Pixel Area in Steradians
        flux_cat('numpy.ndarray'): flux column from catalog
        fluxerr('numpy.ndarray'): corresponded flux error from catalog
        ID('numpy.ndarray'): • ID - (slong), Unique ID of source
        wavelength('float'): Pivot Wavelength in micrometers
        zpvega('float'): vega zeropoint
        mjsr('float'): Flux density (MJy/steradian) producing 1 cps
        pix_sr('float'): pixel area in steradians
        star_col('numpy.ndarray'): classification of stars or galaxies in catalog
        star_or_gal(0 or 1): 0-only stars, 1-only galaxies

    Returns:
        sources_list_cat('list'): list of AB magnitudes of sources
        sources_list_cat_vega('list'): list of vega magnitudes of sources
        sb_list('list'): list of surface brightness of sources in nW m^-2 sr^-1
        err_list('list'): list of error of surface brightness of sources in nW m^-2 sr^-1
    """
    
    #create lists to append magnitudes and surface brightness
    sources_list_cat=[]
    sources_list_cat_vega=[]
    sb_list=[]
    err_list=[]

    #calculate AB magnitude of flux
    flux=(-2.5*np.log10(flux_cat * (10**(-9)) ))+8.9
    #flux vega in nJy
    flux_vega=((mjsr*pix_sr)/ 10**(zpvega/(-2.5)))*(10**15)
    mag_vega=-2.5*np.log10(flux_cat/flux_vega)
    
    
    #over the range of the catalog (all the sources)
    for i in range(len(ID)):
        #calculate surface brightness of source and its error in nW m^-2 sr^-1
        sb=((flux_cat[i]* 10**-9)* ((3*10**8)/(wavelength*10**(-6)))* 10**(-26) * 10**9)/omega_pix
        err=((fluxerr[i]* 10**-9)* ((3*10**8)/(wavelength*10**(-6)))* 10**(-26) * 10**9)/omega_pix
        
        #Append values of AB magnitude, vega magnitude, surface brightness, and error
        #only stars
        if star_or_gal==0:
            
            if star_col[i]==1:
                sources_list_cat.append(flux[i])
                sources_list_cat_vega.append(mag_vega[i])
                sb_list.append(sb)
                err_list.append(err)

        #only galaxies
        if star_or_gal==1:
            
            if star_col[i]!=1:
                sources_list_cat.append(flux[i])
                sources_list_cat_vega.append(mag_vega[i])
                sb_list.append(sb)
                err_list.append(err)

    return sources_list_cat,sources_list_cat_vega, sb_list,err_list


#Create magnitude bins and sum intensity for trilegal simulations and fluxes
def compare_catalog_trilegal(cat_ab_list,cat_vega_list,sb_list,err_list,tri_data,vega_sirius_Jy,wavelength,omega_pix,ab_or_vega,pix_a2):

    """ compare_catalog_trilegal():
        
    Args:
        cat_ab_list('list'): list of AB magnitudes of sources
        cat_vega_list('list'): list of vega magnitudes of sources
        sb_list('list'): list of surface brightness of sources
        err_list('list'): list of error of surface brightness of sources
        tri_data('numpy.ndarray'): array of trilegal simulations
        vega_sirius_Jy('float'): vega mag zeropoint in Jy   
        wavelength('float'): Pivot Wavelength in micrometers
        omega_pix('float'): pixel area in steradians
        ab_or_vega(0 or 1): 0- use AB magnitudes, 1- use Vega magnitudes
        pix_a2('float'): pixel area in arcsecond^2
        
    Returns:
        bins_intensity('numpy.ndarray'): summed intensity per magnitude bin
        bins_range_min('list'): min values of trilegal intensity per magnitude bin
        bins_range_max('list'): max values of trilegal intensity per magnitude bin
        err('numpy.ndarray'): error of surface brightness per magnitude bin
        bins_cat('list'): counts per degree squared for sources per magnitude bin
        bins_range_min_count('list'): min values of counts per degree squared of trilegal simulations per magnitude bin
        bins_range_max_count('list'): max values of counts per degree squared of trilegal simulations per magnitude bin
    """

    #Use sources in either AB magnitudes or vega magnitudes
    if ab_or_vega==0:
        #ab
        cat_list=cat_ab_list
    elif ab_or_vega==1:
        #vega
        cat_list=cat_vega_list

    #find maximum magnitude of sources +1
    num_cat=int(np.floor(np.ma.masked_invalid(cat_list).max()))+1
    
    #find maximum magnitude of trilegal simulations
    max_tri_mag=0
    #for each simulation 
    for sim in range(len(tri_data)):
        #for magnitude in simulation
        for tri_mag in tri_data[sim]:
            
            #calculate AB magnitude
            if ab_or_vega==0:
                #ab
                Fcat = ((10 ** (tri_mag / (-2.5)))* vega_sirius_Jy)
                mag=(-2.5*np.log10(Fcat))+8.9
            #using vega magnitude
            elif ab_or_vega==1:
                #vega
                mag=tri_mag

            #if this magnitude is bigger
            if mag>max_tri_mag:
                #assign max_tri_mag as this magnitude +1
                max_tri_mag=int(np.floor(mag)) +1
    
    #number of bins to use is the biggest magnitude- either from catalog or trilegal
    bins_num= max(num_cat,max_tri_mag)
    
    #create arrays of zeros
    bins=np.zeros(bins_num)
    bins_err=[[] for i in range(len(bins))]
    err=bins.copy()
    bins_intensity=bins.copy()

    #sum intensity per magnitude bin
    for i in range(len(cat_list)):
        #magnitude and error of source
        value=cat_list[i]
        err_value=err_list[i]

        #if magnitude is not nan or inf
        if np.isfinite(value)==True:
            #add count to magnitude bin
            bins[int(np.floor(value))]+=1
            #add surface brightness of source to magnitude bin
            bins_intensity[int(np.floor(value))]+=sb_list[i]
            #append error value to magnitude bin
            bins_err[int(np.floor(value))].insert(int(np.floor(value)),err_value)

    #create list of counts per degree squared per magnitude bin
    bins_cat=[count/pix_a2 for count in bins]

    #for each magnitude bin
    for i in range(len(bins_err)):
        err_sum=0
        #for each error 
        for j in bins_err[i]:
            #if error is not inf or nan
            if np.isfinite(j)==True:
                err_sum+= j**2
        #calculate error of magnitude bin
        calc_err=np.sqrt(err_sum)
        #append calculated error
        err[i]+=calc_err
    
    #create lists for intensity,flux, and ab magnitude for trilegal sims
    intensity_full=[]
    flux_jy_full=[]
    trilegal_full=[]

    #over range of all 10 simulations
    for sim in range(len(tri_data)):
        lIlcat_full=[]
        trilegal_flux_ab=[]
        trilegal_mags=[]
        flux_list=[]

        #For each magnitude in the simulation
        for mag in tri_data[sim]:
            if np.isnan(mag)==False:
                # Flux calculation in Jy
                Fcat = ((10 ** (mag / (-2.5)))* vega_sirius_Jy)
                flux_list.append(Fcat)
                #AB magnitude
                ab=(-2.5*np.log10(Fcat))+8.9

                #append ab magnitude or vega magnitude
                if ab_or_vega==0:
                    trilegal_flux_ab.append(ab)
                elif ab_or_vega==1:
                    trilegal_mags.append(mag)

                #calculate surface brightness in nW m^-2 sr^-1
                lIlcat= ((Fcat)* ((3*10**8)/(wavelength*10**(-6)))* 10**(-26) * 10**9)/omega_pix
                lIlcat_full.append(lIlcat)

        #append list of surface brightness, fluxes, and magnitudes
        intensity_full.append(lIlcat_full)
        flux_jy_full.append(flux_list)
        if ab_or_vega==0:
            trilegal_full.append(trilegal_flux_ab)
        elif ab_or_vega==1:
            trilegal_full.append(trilegal_mags)

    #create bins for trilegal intensity
    full_bins=[]
    full_bins_count=[]

    #over range of all 10 simulations
    for i in range(len(trilegal_full)):
        #create array of zeros
        bins_tri=np.zeros(bins_num)
        bins_intensity_tri=bins_tri.copy()

        #for this simulation
        sim=trilegal_full[i]
        #for each magnitude in this simulation
        for j in range(len(sim)):
            #magnitude
            value=sim[j]
            #if magnitude is not nan
            if np.isfinite(value)==True:
                #rounded magnitude
                index=int(np.floor(value))
                #add count to magnitude bin
                bins_tri[index]+=1
                #add surface brightness to magnitude bin
                bins_intensity_tri[index]+=intensity_full[i][j]
        
        #append array of surface brightness
        full_bins.append(bins_intensity_tri)
        #calculate count per degree squared per magnitude bin
        bins_trilegal=[count/pix_a2 for count in bins_tri]
        full_bins_count.append(bins_trilegal)
    
    #create bins for min and max intensity & count per degree squared for all trilegal simulations
    bins_range_min=[]
    bins_range_max=[]
    bins_range_min_count=[]
    bins_range_max_count=[]

    #per magnitude bin
    for j in range(len(full_bins[0])):
        
        compare_list=[]
        compare_list_count=[]
        #per magnitude bin
        for i in range(len(full_bins)):

            #for one of simulations
            sim=full_bins[i]
            sim_count=full_bins_count[i]
            #value belonging to magnitude bin
            value=sim[j]
            value_count=sim_count[j]
            #assign value of nan if value==0
            if value==0:
                value=np.nan
            if value_count==0:
                value_count=np.nan
            #append value to list
            compare_list.append(value)
            compare_list_count.append(value_count)
        
        #append min and max value of surface brightness and count per degree squared per magnitude bin
        min_val=np.nanmin(compare_list)
        max_val=np.nanmax(compare_list)
        bins_range_min.append(min_val)
        bins_range_max.append(max_val)
        min_count=np.nanmin(compare_list_count)
        max_count=np.nanmax(compare_list_count)
        bins_range_min_count.append(min_count)
        bins_range_max_count.append(max_count)
    
    return bins_intensity,bins_range_min,bins_range_max,err,bins_cat,bins_range_min_count,bins_range_max_count
