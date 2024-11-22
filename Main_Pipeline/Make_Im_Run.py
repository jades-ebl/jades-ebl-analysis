""" Make_Im_Run.py
    
    runs Make_Im.py with variables from MAIN_JADES_Driver.py

    9/19/24
"""

#Imports
#from MAIN_JADES_config import *
from MAIN_JADES_Driver import Read_JADES,ReadJADESCatalog,Read_NircamPSF,Read_TrilegalPickle
#Import function to calculate magnitude zeropoint
from JADES_calc_zeropoint import calc_zeropoint
from JADES_ISL import *
from Make_Im import *
import pickle


def run_make_im(Im_Dir,Cat_Dir,PSF_Dir,Trilegal_Pic_Dir,file_path_zeropoints,Ver,Flux_Col,Flux_Ext,filt_list,pickle_save,flag,flag_ab_or_vega,flag_star_or_gal):
    """ run_make_im():

    Args:
        Im_Dir('string'): full path to filter files of specific version
        Cat_Dir('string'): full path to catalog file of specific version
        PSF_Dir('string'): full path to psf files
        Trilegal_Pic_Dic('string'): full path to trilegal pickle files
        file_path_zeropoints('string'): file path to zeropoint file
        Ver('string'): specifies data version read ('v1' or 'v2' or 'v3')
        Flux_Col('string'): flux column used from catalog ('CIRC0' - Flux of source within circular aperture of 80% enclosed energy radius)
        Flux_Ext(0, 1, or 2): which extension to use flux column (0- circular aperture, 1- background subtracted, 2- psf convolved to f444w)
        filt_list('list'): list of strings specifying which filters to run ('f090w','f115w',...)
        pickle_save('string'): path to save pickle file to
        flag(0 or 1): 0- use masking algorithm applied to noiseless simulation to sum flux of each source, 1- use catalog for the flux of each source
        flag_ab_or_vega(0 or 1): 0- AB magnitude bins, 1- Vega magnitude bins
        flag_star_or_gal(0 or 1): 0- only stars, 1- only galaxies

    Returns:
        cat_bins_filt('numpy.ndarray'): summed intensity of each magnitude bin
        min_bins_filt('list'): minimum intensity of trilegal simulations of each magnitude bin
        max_bins_filt('list'): maximum intensity of trilegal simulations of each magnitude bin
        err_bins_filt('numpy.ndarray'): error of summed intensity of each magnitude bin
        count_cat('numpy.ndarray'): counts per degree squared of each magnitude bin
        count_tri_min('list'): minimum counts per degree squared of trilegal simulations of each magnitude bin
        count_tri_max('list'): maximum counts per degree squared of trilegal simulations of each magnitude bin
        
    """
    save_list=[]

    for Filt in filt_list:
        all_bins=[]
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

        filt_list_full=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f430m','f444w','f460m','f480m']
        wave_list_full=[.901,1.154,1.501,1.845,1.990,2.093,2.786,3.365,3.563,4.092,4.280,4.421,4.624,4.834]

        in_list=filt_list_full.index(Filt)

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

        #Pixel to Arcsecond Ratio
        Pix_Scale= np.sqrt(Pix_A2)

        #Define Frequency (c/λ)
        #nu= ((3*10**8)/(Wave*10**(-6)))

        if Ver=='v2':
            #For these filters there are areas masked out in nosource_mask that have data but no sources in catalog
            if Filt==('f090w' or 'f115w' or 'f150w'):
                with open(r'C:\Users\kasia\nh_pipeline\nosource_mask.pickle','rb') as file:
                    nosource_mask=pickle.load(file)
            else:
                nosource_mask=np.ones(np.shape(im_data))
        else:
            nosource_mask=np.ones(np.shape(im_data))
        
        #sum flux of each source with masking algorithm applied to noiseless simulation
        if flag==0:
            
            #place sources - psf wings
            sourceMap=input_mapper(len(ID_source), im_data.shape,y,x,flux,flux_err,0)[0]
            #convolve sources -psf wings
            sourceMap_convolved= ISL_convolve(psf,sourceMap)
            
            if flag_star_or_gal==0:
                #star
                ab_sources_list,vega_sources_list,sb_list_filt,err_list_filt=Sum_Flux_Mask(Pix_SR,sourceMap_convolved,im_data,err_data,Max_Mag,flux,flux_err,Pix_Scale,ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,Wave,nosource_mask,sirius_vega_zp,mjsr_mean,pix_sr_mean,flag_st,0)
            elif flag_star_or_gal==1:
                #gal
                ab_sources_list,vega_sources_list,sb_list_filt,err_list_filt=Sum_Flux_Mask(Pix_SR,sourceMap_convolved,im_data,err_data,Max_Mag,flux,flux_err,Pix_Scale,ID_source,x,y,s_maj,s_min,theta,bbox_xmin,bbox_xmax,bbox_ymin,bbox_ymax,npix,Wave,nosource_mask,sirius_vega_zp,mjsr_mean,pix_sr_mean,flag_st,1)
        
        #use catalog for flux of each source
        elif flag==1:
            #catalog
            if flag_star_or_gal==0:
                #star
                ab_sources_list,vega_sources_list,sb_list_filt,err_list_filt=Sum_Flux_Mask_cat(Pix_SR,flux,flux_err,ID_source,Wave,sirius_vega_zp,mjsr_mean,pix_sr_mean,flag_st,0)
            elif flag_star_or_gal==1:
                #gal
                ab_sources_list,vega_sources_list,sb_list_filt,err_list_filt=Sum_Flux_Mask_cat(Pix_SR,flux,flux_err,ID_source,Wave,sirius_vega_zp,mjsr_mean,pix_sr_mean,flag_st,1)
        
        #Create magnitude bins and sum intensity for trilegal simulations and fluxes
        cat_bins_filt,min_bins_filt,max_bins_filt,err_bins_filt,count_cat,count_tri_min,count_tri_max=compare_catalog_trilegal(ab_sources_list,vega_sources_list,sb_list_filt,err_list_filt,trilegal_sims,zeropoint,Wave,Pix_SR,flag_ab_or_vega,Pix_A2)
        all_bins.append(Filt)
        all_bins.append(cat_bins_filt)
        all_bins.append(min_bins_filt)
        all_bins.append(max_bins_filt)
        all_bins.append(err_bins_filt)
        all_bins.append(count_cat)
        all_bins.append(count_tri_min)
        all_bins.append(count_tri_max)
        save_list.append(all_bins)
        print('---------------------------------------------------------------------------------------------------------------------------')
    
    with open(pickle_save,'wb') as file:
        pickle.dump(save_list,file)
    return

#For v2------------------
Ver='v2'
flux_col='CIRC0'
flux_ext=0
no_negpix_flag=0
im_dir="D:\JADES\Images2"
#Catalog Directory
Cat_Dir="D:\JADES"
#PSF Directory
PSF_Dir="D:\JADES\PSF"
#Trilegal Pickle Files Directory
Trilegal_Pic_Dir="D:\JADES\ISL\Trilegal_pickles"
file_path_zeropoints = r'D:\JADES\ISL\NRC_ZPs_1126pmap.txt'
filt_list=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f430m','f444w','f460m','f480m']
#-------------------------

#For v1-------------------
# Ver='v1'
# flux_col='CIRC0'
# flux_ext=0
# no_negpix_flag=0
# im_dir="D:\JADES\Images"
# #Catalog Directory
# Cat_Dir="D:\JADES"
# #PSF Directory
# PSF_Dir="D:\JADES\PSF"
# #Trilegal Pickle Files Directory
# Trilegal_Pic_Dir="D:\JADES\ISL\Trilegal_pickles"
# file_path_zeropoints = r'D:\JADES\ISL\NRC_ZPs_1126pmap.txt'
# filt_list=['f090w','f115w','f150w','f200w','f277w','f335m','f356w','f410m','f444w']
#---------------------------

#For v3------------------
# Ver='v3'
# flux_col='CIRC0'
# flux_ext=0
# no_negpix_flag=0
# im_dir="D:\JADES\Images3"
# #Catalog Directory
# Cat_Dir="D:\JADES"
# #PSF Directory
# PSF_Dir="D:\JADES\PSF"
# #Trilegal Pickle Files Directory
# Trilegal_Pic_Dir="D:\JADES\ISL\Trilegal_pickles"
# file_path_zeropoints = r'D:\JADES\ISL\NRC_ZPs_1126pmap.txt'
# filt_list=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f444w']
#-------------------------

#save path
pickle_save=r'd:\JADES\ISL\Trilegal_pickles\v2_ab_gals_makeim.pickle'
#0-masking,1-catalog
flag=1
#0-ab,1-vega
flag_ab_or_vega=0
#0-stars,1-gal
flag_star_or_gal=1

run_make_im(im_dir,Cat_Dir,PSF_Dir,Trilegal_Pic_Dir,file_path_zeropoints,Ver,flux_col,flux_ext,filt_list,pickle_save,flag,flag_ab_or_vega,flag_star_or_gal)