
""" JADES_config.py

    Sets variables for JADES_driver

    2/19/24
"""

#Image Directory for specific version
Im_Dir="D:\JADES\Images"
#Catalog Directory
Cat_Dir="D:\JADES"

""" Filters for each Data version
v1:
    F090W, F115W, F150W,        F200W,        F277W, F335M, F356W, F410M,        F444W
v2:
    F090W, F115W, F150W, F182M, F200W, F210M, F277W, F335M, F356W, F410M, F430M, F444W, F460M, F480M
"""
#Filter (lowercase)
Filt='f090w'
#Version (lowercase: 'v1' or 'v2')
Ver='v1'


#AB Mag Column name from Catalog (Ex: 'CIRC0', uppercase)
""" 
    CIRC0- Flux of source within circular aperture of 80% enclosed energy radius
    CIRC1- Flux of source within circular aperture of 0.10" radius
    CIRC2- Flux of source within circular aperture of 0.15" radius
    CIRC3- Flux of source within circular aperture of 0.25" radius
    CIRC4- Flux of source within circular aperture of 0.30" radius
    CIRC5- Flux of source within circular aperture of 0.35" radius
    CIRC6- Flux of source within circular aperture of 0.50" radius
"""
Flux = 'CIRC0'
#AB Mag Table Ext: T4: Circ aperture (0), T5: Bsub (Background-Subtracted) (1), T6: Conv (psf-convolved to f444w) (2), Default 0
Flux_Col=0
#Flag (0- nothing, 1- Zscale with LinearStretch plot of filter image data, 2- Plot of Mask, 3- Plot of Mask*Data, 4- Histogram, 5-Chunking Plot of Data, 6- Chunking Plot of Mask*Data, 7 Chunk Histogram, 8 ISL)
Flag=0

#Do you want to save figure? 'Y' yes or 'N' no
save_fig= 'N'
#where to save figure (path)
file_save = 'd:\JADES\Chunking\Chunk_1_hist.pdf'
#dpi setting for figure
dpi_num=300

#Bounding Limits for Chunking Box
y_min=6037
y_max=7622
x_min=7448
x_max=9197
