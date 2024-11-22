# Before running this code you MUST edit /home/symons/anaconda3/lib/python3.7/site-packages/astrobase/services/trilegal.py
# ADD to TRILEGAL_FILTER_SYSTEMS dict:
#     'gaiaDR2': {
#         'desc': "Gaia's DR2 G, G_BP, G_RP (Vegamags, Gaia passbands from Evans et al. 2018)",
#         'table': 'tab_mag_odfnew/tab_mag_gaiaDR2.dat'
#     },
# COMMENT OUT extinction_info and Av_infinity lines w/ try/except ect. in query_galcoords function and replace them with:
#     #instead, use website default
#     Av_infinity = 0.0378
#     inputparams['extinction_infty'] = '%.5f' % Av_infinity

import os
import gzip
import shutil
import astrobase.services.trilegal
from JADES_getcoords import JADES_getcoords


dir_base = r'd:\JADES\ISL\trilegal_v3' #base directory for outputs

#For v1-------------------
#ver='v1'
#im_dir="D:\JADES\Images"
#filt_list_full=['F090W','F115W','F150W','F200W','F277W','F335M','F356W','F410M','F444W']
#---------------------------
#For v2------------------
#ver='v2'
#im_dir="D:\JADES\Images2"
#filt_list_full=['F090W', 'F115W', 'F150W', 'F182M', 'F200W', 'F210M', 'F277W', 'F335M', 'F356W', 'F410M', 'F430M', 'F444W', 'F460M', 'F480M']
#-------------------------
#For v3------------------
ver='v3'
im_dir="D:\JADES\Images3"
filt_list_full=['F090W','F115W','F150W','F182M','F200W','F210M','F277W','F335M','F356W','F410M','F444W']
#-------------------------

filters,field_RA,field_DEC,field_degs=JADES_getcoords(im_dir,ver,filt_list_full)

field_fieldNums = [i for i in range(1,len(field_DEC)+1)] #field num array, make a list if not going 1 to lenDEC
field_numPerField = 10 #num times to run the same field

field_numFields = len(field_DEC) #get len now
if( (len(field_DEC) != len(field_RA)) | (len(field_DEC) != len(field_fieldNums)) | (len(field_degs) != len(field_DEC)) ):
    print('ERROR: Len field_DEC is '+str(len(field_DEC))+', len field_RA is '+str(len(field_RA))+', len field_fieldNum is '+str(len(field_fieldNums))+', len field_degs is '+str(len(field_degs))+' and at least one of them doesn\'t match.')
    import sys
    sys.crash() #not a real call

#END IF
    
for i in range(0,field_numFields):
    dir_curr = os.path.join(dir_base,str(filters[i]).zfill(2)) #create current path
    if os.path.isdir(dir_curr) == False:
        os.makedirs(dir_curr) #create dir if not a dir
    #END IF

    for j in range(0,field_numPerField):
        return_dict = astrobase.services.trilegal.query_radecl(field_RA[i], field_DEC[i], filtersystem='jwst_nircam', field_deg2= field_degs[i],
                    usebinaries=True, extinction_sigma=0, magnitude_limit=32.0, maglim_filtercol=1, trilegal_version=1.6,
                    extraparams=None, forcefetch=True, cachedir=dir_curr, verbose=True,
                    timeout=60.0, refresh=60.0, maxtimeout=1500.0) #read trilegal file, it is zipped


        with gzip.open(return_dict['tablefile'], 'rb') as f_in:
            with open(return_dict['tablefile'][:return_dict['tablefile'].find('.')]+str(j)+'.dat', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out) #un gzips and saves file, inspired by https://stackoverflow.com/a/44712152
            #END WITH
        #END WITH
    #END FOR j
    os.remove(return_dict['tablefile']) #cleanup (.gz file name is same for all FOR j runs, so only 1 .gz to cleanup at the end)
#END FOR i

print('done')
