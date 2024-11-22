""" JADES_runDriver.py
    #Runs pipeline for all filters of a specified data version, and saves measurement numbers to csv and pickle file
    11/22/24
"""

#Imports
from astropy.io import fits
from astropy.io import ascii
import os
from astropy.table import Table
import astropy.wcs as wcs
import numpy as np
#import Masking Function
from JADES_Mask import *
#Import PSF ISL Function
from JADES_ISL import *
#Import function to calculate magnitude zeropoint
from JADES_calc_zeropoint import *
#Import faint ISL Function
from JADES_trilegalisl import *


def run_Driver_allfilt(Im_Dir, Cat_Dir, PSF_Dir, Trilegal_Pic_Dir, file_path_zeropoints, Ver, Flux_Col, Flux_Ext, this_filt_list,comment_input,save_path,pickle_save,no_negpix_flag):
    """run_Driver_allfilt():
        #Runs pipeline for all filters of a specified data version, and saves measurement numbers to csv and pickle file
    Args:
        Im_Dir ('string'): full path to filter files of specific version
        Cat_Dir ('string'): full path to catalog file of specific version
        PSF_Dir ('string'): full path to psf files
        Trilegal_Pic_Dir ('string'): full path to trilegal pickle files
        file_path_zeropoints ('string'): file path to zeropoint file
        Ver ('string'): specifies data version read ('v1' or 'v2' or 'v3')
        Flux_Col ('string'): flux column used from catalog ('CIRC0' - Flux of source within circular aperture of 80% enclosed energy radius)
        Flux_Ext (0, 1, or 2): which extension to use flux column (0- circular aperture, 1- background subtracted, 2- psf convolved to f444w)
        this_filt_list ('list'): list of strings specifying which filters to run ('f090w','f115w',...)
        comment_input ('list'): list of strings to add comments to saved file
        save_path ('string'): path to save csv file to
        pickle_save ('string'): path to save pickle file to
        no_negpix_flag (0 or 1): 1- remove negative pixels in data in masking function measurement, 0- don't remove negative pixels
    """
    print('---------------------------------------------------------------------------------------------------------------------------')
    from MAIN_JADES_Driver import Read_JADES, ReadJADESCatalog, Read_NircamPSF, Read_TrilegalPickle
    #set up table column names and data types
    save_table=Table(names=('Filter','Wavelength','Percent_Pix','Percent_MaskErr','Percent_ScaleMag','Percent_ScaleArea','I*','err_I*','Isky','err_Isky','Igal','err_Igal','Ifaint','err_Ifaint','IPSF','err_IPSF','err_Isky(phot_cal)','ICIB','err_ICIB'),dtype=(str,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))
    save_table.meta['comments'] = comment_input

    #For each filter
    for Filt in this_filt_list:
        
        #full list of filters and wavelengths for JADES
        filt_list_full=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f430m','f444w','f460m','f480m']
        wave_list_full=[.901,1.154,1.501,1.845,1.990,2.093,2.786,3.365,3.563,4.092,4.280,4.421,4.624,4.834]

        in_list=filt_list_full.index(Filt)
        #Match Pivot Wavelength to Filter
        Wave=wave_list_full[in_list]
        

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

        #Pixel to Arcsecond Ratio
        Pix_Scale= np.sqrt(Pix_A2)

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

        #count percent of sources to go through each scaling method & how many sources have masked error
        count_maskerr=0
        count_scalemag=0
        count_scalearea=0
        len_sources=0
        for i in range(len(ID_source)):
            #If source in catalog is in valid data area
            if err_data[int(np.round(y[i])),int(np.round(x[i]))]==True:
                #count number of sources included within data filter
                len_sources+=1
                #count sources with masked error
                if ma.is_masked(flux_err[i])==True:
                    count_maskerr+=1
                #calculate AB magnitude of flux
                thismag=(-2.5*np.log10(flux[i] * (10**(-9)) ))+8.9

                #if the magnitude is a valid number (not infinite or nan) - magnitude scaling
                if np.abs(thismag) != np.inf and np.isnan(thismag)!=True:
                    count_scalemag+=1

                #if magnitude is inf or nan - area scaling     
                else:
                    count_scalearea+=1
        print('All sources=',len(ID_source))
        print('# Sources=',len_sources)
        percent_maskerr=(count_maskerr/len_sources)*100
        percent_scalemag=(count_scalemag/len_sources)*100
        percent_scalearea=(count_scalearea/len_sources)*100
        print('Percent of sources with masked error =',percent_maskerr)
        print('Percent of sources for first scaling method (mag) =',percent_scalemag)
        print('Percent of sources for second scaling method (area) =',percent_scalearea)

        print('---------------------------------------------------------------------------------------------------------------------------')

        #Calculate number of pixels
        data_area=im_data*nosource_mask
        num_pixels=np.nansum(np.isfinite(np.where(err_data,data_area,np.nan)))

        #set trilegal mag mask limit to 27 ab mag
        max_mask=27

        #Mask Image Data 
        mask,percent_pix,hist,mean_nonmask,mean_mask,mean_mask_gal_plus_dimstars=Final_Masking_new(im_data,err_data,Max_Mag,flux,flux_err,Pix_Scale,ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,0,Wave,nosource_mask,flag_st,no_negpix_flag,max_mask)
        mask_err,err_percent_pix,hist,mean_nonmask_err,mean_mask_err,mean_mask_gal_plus_dimstars_err=Final_Masking_new(im_data,err_data,Max_Mag,flux,flux_err,Pix_Scale,ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,1,Wave,nosource_mask,flag_st,no_negpix_flag,max_mask)
        print('---------------------------------------------------------------------------------------------------------------------------')
        
        print("MEAN mask- all")
        print("Mean surface brightness of resolved stars and galaxies = ", mean_mask, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Mean surface brightness after masking stars and galaxies= ", mean_nonmask, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        err_resolved=np.abs(mean_mask_err-mean_mask)
        err_nonmask= np.abs(mean_nonmask_err-mean_nonmask)
        print("Error of surface brightness of resolved stars and galaxies = ", err_resolved, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        print("Error of surface brightness after masking stars and galaxies= ", err_nonmask, " nW m^-2 sr^-1", "(", Filt, Ver, ")")
        
        
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
        
        #save measurements to table
        save_table.add_row([Filt,Wave,percent_pix,percent_maskerr,percent_scalemag,percent_scalearea,mean_mask,err_resolved,mean_nonmask,err_nonmask,mean_mask_gal_plus_dimstars,err_resolved_gal_plus_dimstars,mean_tri,err_tri,ISL_mean,err_psf,err_photcal,lam_CIB,Quad_Error])
        print('---------------------------------------------------------------------------------------------------------------------------')
    print(save_table)

    #save table to csv file
    ascii.write(save_table,str(save_path), overwrite=True)
    #save table to pickle file
    with open(pickle_save,'wb') as file:
        pickle.dump(save_table,file) 


#For v1-------------------
# this_ver='v1'
# flux_col='CIRC0'
# flux_ext=0
# no_negpix_flag=0
# #Filter Data Directory
# im_dir="D:\JADES\Images"
# #Catalog Directory
# cat_dir="D:\JADES"
# #PSF Directory
# psf_dir="D:\JADES\PSF"
# #Trilegal Pickle Files Directory
# Trilegal_Pic_Dir="D:\JADES\ISL\Trilegal_pickles"
# #Zeropoint file path
# file_path_zeropoints = r'D:\JADES\ISL\NRC_ZPs_1126pmap.txt'
# #v1 filters
# filt_list_full=['f090w','f115w','f150w','f200w','f277w','f335m','f356w','f410m','f444w']
#---------------------------

#For v2------------------
# this_ver='v2'
# flux_col='CIRC0'
# flux_ext=0
# no_negpix_flag=0
# #Filter Data Directory
# im_dir="D:\JADES\Images2"
# #Catalog Directory
# cat_dir="D:\JADES"
# #PSF Directory
# psf_dir="D:\JADES\PSF"
# #Trilegal Pickle Files Directory
# Trilegal_Pic_Dir="D:\JADES\ISL\Trilegal_pickles"
# #Zeropoint file path
# file_path_zeropoints = r'D:\JADES\ISL\NRC_ZPs_1126pmap.txt'
# #v2 filters
# filt_list_full=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f430m','f444w','f460m','f480m']
#-------------------------

#For v3------------------
this_ver='v3'
flux_col='CIRC0'
flux_ext=0
no_negpix_flag=0
#Filter Data Directory
im_dir="D:\JADES\Images3"
#Catalog Directory
cat_dir="D:\JADES"
#PSF Directory
psf_dir="D:\JADES\PSF"
#Trilegal Pickle Files Directory
Trilegal_Pic_Dir="D:\JADES\ISL\Trilegal_pickles"
#Zeropoint file path
file_path_zeropoints = r'D:\JADES\ISL\NRC_ZPs_1126pmap.txt'
#v3 filters
filt_list_full=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f444w']
#-------------------------

#comments and save path
this_comment_input=['Ver 3 (11/21/24)','using: circular aperture ext'] 
save_path_csv=r'D:\JADES\Calculations_New\v3_circap_1.csv'
save_pickle_path=r'D:\JADES\Calculations_New\v3_circap_1.pickle'

run_Driver_allfilt(im_dir, cat_dir, psf_dir, Trilegal_Pic_Dir, file_path_zeropoints, this_ver, flux_col, flux_ext, filt_list_full,this_comment_input,save_path_csv,save_pickle_path,no_negpix_flag)
