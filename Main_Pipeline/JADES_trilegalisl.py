
""" JADES_trilegalisl.py

    Calculate surface brightness of stars over masking limit with trilegal simulations

    10/7/24
"""

#Imports
import numpy as np

def jades_trilegalisl(data, num_pixels,PIXAR_A2, vega_sirius_Jy, nu, max_mask):
    """jades_trilegalisl(data, num_pixels,PIXAR_A2, vega_sirius_Jy, nu, max_mask)

    Args:
        data (numpy.ndarray): trilegal simulations
        num_pixels ('float'): number of pixels of JADES image
        PIXAR_A2 ('float'): Nominal pixel area in arcsec^2
        vega_sirius_Jy ('float'): zeropoint vega_sirius_Jy
        nu ('float'): frequency (c/λ)
        max_mask ('float'): limit of masking stars in AB Mag

    Returns:
        trilegal_mean ('float'): Surface brightness of stars above masking limit  
        err_full: _description_
    """

    #compute the area of a JADES image
    surveyarea = num_pixels * (PIXAR_A2/ 3600**2)

    #number of simulations (10)
    nfiles = np.shape(data)[0]

    #create list to append lIlcat
    lIlcat_list=[]

    #For each simulation
    for jfile in range(nfiles):
        sim = data[jfile]
        #create list to append Fcat
        flux_list=[]

        #For each magnitude in the simulation
        for mag in sim:
            if np.isnan(mag)==False:
                #convert from mag to flux
                #Flux calculation in Jy
                Fcat = ((10 ** (mag / (-2.5)))* vega_sirius_Jy)
                
                #AB magnitude
                ab=(-2.5*np.log10(Fcat))+8.9
                #if valid magnitude
                if np.isfinite(ab)==True:
                    #if magnitude is above masking limit append flux
                    if ab>max_mask:
                        flux_list.append(Fcat)

        #calculate surface brightness of stars in nW m^-2 sr^-1
        lIlcat = (10**(-26) * 10**(9) * nu * np.sum(flux_list) )/ (surveyarea * (np.pi / 180) ** 2) # this one is ISL - per radian squared
        lIlcat_list.append(lIlcat)
        
    #calculate mean and error of simulations
    trilegal_mean=np.mean(lIlcat_list)
    trilegal_err=np.std(lIlcat_list)

    return trilegal_mean, trilegal_err

def jades_trilegalisl_chunk(data, num_pixels, num_pix_chunk, PIXAR_A2, vega_sirius_Jy, nu, max_mask, chunk_num):
    """jades_trilegalisl_chunk(data, num_pixels, num_pix_chunk, PIXAR_A2, vega_sirius_Jy, nu, max_mask, chunk_num)

    Args:
        data (numpy.ndarray): trilegal simulations
        num_pixels ('float'): number of pixels of JADES image
        num_pix_chunk ('float'): number of pixels in chunk
        PIXAR_A2 ('float'): Nominal pixel area in arcsec^2
        vega_sirius_Jy ('float'): zeropoint vega_sirius_Jy
        nu ('float'): frequency (c/λ)
        max_mask ('float'): limit of masking stars in AB Mag
        chunk_num ('int'): chunk number

    Returns:
        trilegal_mean ('float'): Surface brightness of stars above masking limit  
        err_full: _description_
    """

    #compute the area of a JADES image
    surveyarea = num_pixels * (PIXAR_A2/ 3600**2)
    #compute percent of pixels in chunk to image
    percent_area= num_pix_chunk/num_pixels

    #number of simulations (10)
    nfiles = np.shape(data)[0]

    #create list to append lIlcat
    lIlcat_list=[]

    #For each simulation
    for jfile in range(nfiles):
        sim = data[jfile]
        #number of stars simulated 
        len_sim= len(np.where(np.isnan(sim))[0])
        #number of stars in chunk
        num_stars= np.round(len_sim* percent_area)
        
        #return random indicies simulated in chunk
        rng=np.random.default_rng(seed=chunk_num)
        rints=rng.integers(low=0,high=len_sim,size=int(num_stars))
        #create list to append Fcat
        flux_list=[]

        #For each random integer
        for ind in rints:
            mag = sim[ind]
            if np.isnan(mag)==False:
                #convert from mag to flux
                #Flux calculation in Jy
                Fcat = ((10 ** (mag / (-2.5)))* vega_sirius_Jy)
                
                #AB magnitude
                ab=(-2.5*np.log10(Fcat))+8.9
                #if valid magnitude
                if np.isfinite(ab)==True:
                    #if magnitude is above masking limit append flux
                    if ab>max_mask:
                        flux_list.append(Fcat)

        #calculate surface brightness of stars in nW m^-2 sr^-1
        lIlcat = (10**(-26) * 10**(9) * nu * np.sum(flux_list) )/ (surveyarea * (np.pi / 180) ** 2) # this one is ISL - per radian squared
        lIlcat_list.append(lIlcat)
        
    #calculate mean and error of simulations
    trilegal_mean=np.mean(lIlcat_list)
    trilegal_err=np.std(lIlcat_list)

    return trilegal_mean, trilegal_err

