"""JADES_getcoords.py

    get filters list, field_DEC, field_RA, and field_degs for filter version

    10/11/24
"""
#Imports
import pickle
import numpy as np
from MAIN_JADES_Driver import Read_JADES

def JADES_getcoords(Im_Dir,Ver,filt_list):
    """JADES_getcoords()

    Args:
        Im_Dir ('string'): full path to filter files of specific version
        Ver ('string'): specifies data version read ('v1' or 'v2' or 'v3')
        filt_list ('list'): list of strings specifying which filters to run ('f090w','f115w',...)

    Returns:
        filters ('list'):  list of strings specifying which filters ran
        field_RA ('list'): list of targ RA of each filter
        field_DEC ('list'): list of targ DEC of each filter
        field_degs ('list'): list of survey area in square degrees of each filter
    """

    #create lists to append values
    filters=[]
    field_RA=[]
    field_DEC=[]
    field_degs=[]

    #For each filter
    for Filt in filt_list:
        #Read in filter Data
        im_data,err_data,photmjsr,Pix_SR,Pix_A2,Targ_RA,Targ_DEC=Read_JADES(Im_Dir,Filt,Ver)
        print('---------------------------------------------------------------------------------------------------------------------------')
        #Multiply data by nosource mask
        if Ver=='v2':
            #For these filters there are areas masked out in nosource_mask that have data but no sources in catalog
            if Filt==('f090w' or 'f115w' or 'f150w'):
                with open(r'C:\Users\kasia\nh_pipeline\nosource_mask.pickle','rb') as file:
                    nosource_mask=pickle.load(file)
            else:
                nosource_mask=np.ones(np.shape(im_data))
        else:
            nosource_mask=np.ones(np.shape(im_data))

        #calculate survey area in square degrees
        data_area=im_data * nosource_mask
        num_pixels=np.sum(np.isfinite(np.where(err_data,data_area,np.nan)))
        survey_area= num_pixels * (Pix_A2/3600**2)

        #append values
        filters.append(Filt)
        field_RA.append(Targ_RA)
        field_DEC.append(Targ_DEC)
        field_degs.append(survey_area)
    
    
    return filters,field_RA,field_DEC,field_degs


#For v2------------------
#ver='v2'
#im_dir="D:\JADES\Images2"
#filt_list_full=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f430m','f444w','f460m','f480m']
#-------------------------

#For v1-------------------
#ver='v1'
#im_dir="D:\JADES\Images"
#filt_list_full=['f090w','f115w','f150w','f200w','f277w','f335m','f356w','f410m','f444w']
#---------------------------

#For v3------------------
#ver='v3'
#im_dir="D:\JADES\Images3"
#filt_list_full=['f090w','f115w','f150w','f182m','f200w','f210m','f277w','f335m','f356w','f410m','f444w']
#-------------------------

#filts,ras,decs,degs=JADES_getcoords(im_dir,ver,filt_list_full)