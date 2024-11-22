""" JADES_calc_zeropoint.py

    Calculates zero point of magnitude

    6/8/24
"""
import numpy as np


def calc_zeropoint(filt, file_path):
    """calc_zeropoint()

    Args:
        filt ('string'): uppercase string that specifies filter
        file_path ('string'): file path to zeropoint file

    Returns:
        zeropoint ('float'): vega-sirius_Jy
        vega_zeropoint ('float'): zp_vega
        sirius_vega_zeropoint ('float'): zp_vega-sirius
        mjsr_mean ('float'): PHOTMJSR
        pix_sr_mean ('float'): mean_pix_sr
    """
    
    #open file
    file=open(file_path)
    print('READING ZEROPOINT FILE')
    #list to append looped file data
    zeropoint_list=[]

    #other lists
    vega_zp_list=[]
    vega_sirius_zeropoints_list=[]
    mjsr_values_list=[]
    mean_pix_sr_list=[]

    #while there are lines to read
    while (True):
        #read each line
        line=file.readline()
            
        #if first line in file
        while line.startswith('#'):
            #read next line
            line=file.readline()

        #break loop at end
        if not line:
            file.close()
            break
            
        #split each line to get values from inputted filter
        filt_name= line.split()[1]
        #zp_vega
        vega_zp=line.split()[7]
        #vega-sirius_Jy
        zeropoints = line.split()[13]
        #zp_vega-sirius
        vega_sirius_zeropoints = line.split()[11]
        #PHOTMJSR
        mjsr_values = line.split()[5]
        #mean_pix_sr
        mean_pix_sr = line.split()[17]

        #append values with specified filter
        if filt_name == str('CLEAR+'+filt):
            zeropoint_list.append(float(zeropoints))
            vega_zp_list.append(float(vega_zp))
            vega_sirius_zeropoints_list.append(float(vega_sirius_zeropoints))
            mjsr_values_list.append(float(mjsr_values))
            mean_pix_sr_list.append(float(mean_pix_sr))
    
    #take mean of values
    zeropoint=np.mean(zeropoint_list)
    vega_zeropoint=np.mean(vega_zp_list)
    sirius_vega_zeropoint=np.mean(vega_sirius_zeropoints_list)
    mjsr_mean=np.mean(mjsr_values_list)
    pix_sr_mean=np.mean(mean_pix_sr_list)
    
    return zeropoint,vega_zeropoint, sirius_vega_zeropoint, mjsr_mean, pix_sr_mean


# filt='F090W'
# file_path = r'd:\JADES\ISL\NRC_ZPs_1126pmap.txt'

#calc_zeropoint(filt,file_path)













