
""" MAIN_JADES_config.py

    Sets variables for MAIN_JADES_driver

    5/12/24
"""


#Image Directory for specific version
Im_Dir="D:\JADES\Images2"
#Catalog Directory
Cat_Dir="D:\JADES"
#PSF Directory
PSF_Dir="D:\JADES\PSF"
#Trilegal Pickle Files Directory
Trilegal_Pic_Dir="D:\JADES\ISL\Trilegal_pickles"
#File path for zeropoints
file_path_zeropoints = r'D:\JADES\ISL\NRC_ZPs_1126pmap.txt'


""" Filters for each Data version
v1:
    F090W, F115W, F150W,        F200W,        F277W, F335M, F356W, F410M,        F444W
v2:
    F090W, F115W, F150W, F182M, F200W, F210M, F277W, F335M, F356W, F410M, F430M, F444W, F460M, F480M

v3 (GOODS_N):
    F090W, F115W, F150W, F182M, F200W, F210M, F277W, F335M, F356W, F410M,      , F444W
"""
#Filter (lowercase)
Filt='f090w'
#Version (lowercase: 'v1' or 'v2' or 'v3')
Ver='v2'

#Flag (1-EBL, 2- EBL for chunk)
Flag=1

#Flux Column name from Catalog (Ex: 'CIRC0', uppercase)
""" 
    CIRC0- Flux of source within circular aperture of 80% enclosed energy radius
    CIRC1- Flux of source within circular aperture of 0.10" radius
    CIRC2- Flux of source within circular aperture of 0.15" radius
    CIRC3- Flux of source within circular aperture of 0.25" radius
    CIRC4- Flux of source within circular aperture of 0.30" radius
    CIRC5- Flux of source within circular aperture of 0.35" radius
    CIRC6- Flux of source within circular aperture of 0.50" radius
"""
Flux_Col = 'CIRC0'
#Flux Ext: T4: Circ aperture (0), T5: Bsub (Background-Subtracted) (1), T6: Conv (psf-convolved to f444w) (2), Default 0
Flux_Ext=0
