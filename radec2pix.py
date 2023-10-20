import numpy as np
import argparse

def radec2pix(alpha, delta, astrometry, varargin):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('alpha', type=float)
    parser.add_argument('delta', type=float)
    parser.add_argument('astrometry', type=dict)
    parser.add_argument('--cdshift', type=float, default=1) # optional argument
    
    # parse the arguments
    args = parser.parse_args()
    
    alpha = args.alpha
    delta = args.delta
    astrometry = args.astrometry
    cdshift = args.cdshift
    
    # set up some basic info to pass to the coordinate rotator
    wcsparams = {
        'lonpole': astrometry['lonpole'],
        'latpole': astrometry['latpole'],
        'crval': [astrometry['crval1'], astrometry['crval2']]
    }
    
    # rotate the coordinates
    xsi, eta = wcssph2xy(alpha, delta, wcsparams)

    # make the plate scale matrix
    cd = cdshift * np.array([[astrometry['cd1_1'], astrometry['cd1_2']],
                         [astrometry['cd2_1'], astrometry['cd2_2']]])

    # calculate the inverse to go to the pixel number
    cdinv = np.linalg.inv(cd)
    
    # multiply by the relevant vectors
    u = cdinv[0, 0] * xsi + cdinv[0, 1] * eta
    v = cdinv[1, 0] * xsi + cdinv[1, 1] * eta
    
    # undo the distortion
    up = u
    vp = v
    
    if 'ap_order' in astrometry:
        
        if astrometry['ap_order'] == 1:
            if 'ap_0_0' in astrometry:
                up = up + astrometry['ap_0_0'] * u ** 0 * v ** 0
            if 'ap_1_0' in astrometry:
                up = up + astrometry['ap_1_0'] * u ** 1 * v ** 0
            if 'ap_0_1' in astrometry:
                up = up + astrometry['ap_0_1'] * u ** 0 * v ** 1
            if 'ap_1_1' in astrometry:
                up = up + astrometry['ap_1_1'] * u ** 1 * v ** 1
                
            
            if 'bp_0_0' in astrometry:
                vp = vp + astrometry['bp_0_0'] * u ** 0 * v ** 0
            if 'bp_1_0' in astrometry:
                vp = vp + astrometry['bp_1_0'] * u ** 1 * v ** 0
            if 'bp_0_1' in astrometry:
                vp = vp + astrometry['bp_0_1'] * u ** 0 * v ** 1
            if 'bp_1_1' in astrometry:
                vp = vp + astrometry['bp_1_1'] * u ** 1 * v ** 1

        elif astrometry['ap_order'] == 2:
            if 'ap_0_0' in astrometry:
                up = up + astrometry['ap_0_0'] * u ** 0 * v ** 0
            if 'ap_1_0' in astrometry:
                up = up + astrometry['ap_1_0'] * u ** 1 * v ** 0
            if 'ap_0_1' in astrometry:
                up = up + astrometry['ap_0_1'] * u ** 0 * v ** 1
            if 'ap_1_1' in astrometry:
                up = up + astrometry['ap_1_1'] * u ** 1 * v ** 1
            if 'ap_0_2' in astrometry:
                up = up + astrometry['ap_0_2'] * u ** 0 * v ** 2
            if 'ap_1_2' in astrometry:
                up = up + astrometry['ap_1_2'] * u ** 1 * v ** 2
            if 'ap_2_0' in astrometry:
                up = up + astrometry['ap_2_0'] * u ** 2 * v ** 0
            if 'ap_2_1' in astrometry:
                up = up + astrometry['ap_2_1'] * u ** 2 * v ** 1
            if 'ap_2_2' in astrometry:
                up = up + astrometry['ap_2_2'] * u ** 2 * v ** 2
                
            if 'bp_0_0' in astrometry:
                vp = vp + astrometry['bp_0_0'] * u ** 0 * v ** 0
            if 'bp_1_0' in astrometry:
                vp = vp + astrometry['bp_1_0'] * u ** 1 * v ** 0
            if 'bp_0_1' in astrometry:
                vp = vp + astrometry['bp_0_1'] * u ** 0 * v ** 1
            if 'bp_1_1' in astrometry:
                vp = vp + astrometry['bp_1_1'] * u ** 1 * v ** 1
            if 'bp_0_2' in astrometry:
                vp = vp + astrometry['bp_0_2'] * u ** 0 * v ** 2
            if 'bp_1_2' in astrometry:
                vp = vp + astrometry['bp_1_2'] * u ** 1 * v ** 2
            if 'bp_2_0' in astrometry:
                vp = vp + astrometry['bp_2_0'] * u ** 2 * v ** 0
            if 'bp_2_1' in astrometry:
                vp = vp + astrometry['bp_2_1'] * u ** 2 * v ** 1
            if 'bp_2_2' in astrometry:
                vp = vp + astrometry['bp_2_2'] * u ** 2 * v ** 2
            
        else:
            iaX, jaY = np.meshgrid(np.arange(astrometry['ap_order'] + 1), np.arange(astrometry['ap_order'] + 1))
            #??
            thisord_mat = [['ap_' + str(jaY[i, j]) + '_' + str(iaX[i, j]) for j in range(jaY.shape[1])] for i in range(jaY.shape[0])]
            for ia in range(astrometry['ap_order'] + 1):
                for ja in range(astrometry['ap_order'] + 1):
                    thisord = thisord_mat[ia][ja]
                    if thisord in astrometry:
                        up = up + astrometry[thisord] * u ** ia * v ** ja
                        
                        
            iaX, jaY = np.meshgrid(np.arange(astrometry['bp_order'] + 1), np.arange(astrometry['bp_order'] + 1))
            #??
            thisord_mat = [['bp_' + str(jaY[i, j]) + '_' + str(iaX[i, j]) for j in range(jaY.shape[1])] for i in range(jaY.shape[0])]
            for ia in range(astrometry['bp_order'] + 1):
                for ja in range(astrometry['bp_order'] + 1):
                    thisord = thisord_mat[ia][ja]
                    if thisord in astrometry:
                        vp = vp + astrometry[thisord] * u ** ia * v ** ja
    
    # add back the reference pixel
    x = vp + astrometry['crpix2']
    y = up + astrometry['crpix1']
    
    return x, y
