import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from skimage import io
import os, fnmatch
from os import listdir
from astropy.io import fits

import reduction_tools
from reduction_tools import tof2l
from reduction_tools import find_nearest
import phasefitting_1D

import time
   
def phase_ratio_linearcomb_2D(Ttof,spectrum,phase1lac,phase2lac,phase_spectrum,lambda_range_norm,lambda_range_edges,filemask=0,auto_mask=True,mask_thresh=[0.05, 0.95],est_phi=0.5,method='least_squares',bool_SG=False,SG_w=5,SG_n=1,bool_save=False,bool_print=False,debug_flag=False,debug_idx=[160,200]):
    """ Performs edge fitting with gaussian model to stack of TOF images (x,y,lambda)
    !!! IMPORTANT the data must be transmission or it won't work due to boundary conditions !!!
    
    INPUTS:
    Ttof = 3d matrix with the stack of TRANSMISSION (I/I0) TOF images (x,y,lambda) 
    spectrum_l = spectrum, length of this ndarray must correspond to size of Ttof(lambda)
    lambda_range = range corresponding to lambda where to perform the fitting
    filemask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_pos = estimated bragg edge position (in spectrum_l dimension)
    est_wid = estimated bragg edge width (in spectrum_l dimension)
    est_h = estimated bragg edge height (in spectrum_l dimension)
    bool_smooth = set to True to perform Savitzky-Golay filtering of the transmission derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_flag = set to True to test on single pixel (will actually also perform the whole image fitting if not stopped)
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    'edge_position' : edge position 
    'edge_height': edge height 
    'edge_width': edge width  
    'edge_slope': edge slope 
    'median_image': median Transmission image in the selected lambda range
    """ 

    if(filemask.any()):
        mymask = filemask
        plt.figure()
        plt.subplot(1,2,1), plt.imshow(np.median(Ttof,axis=2)), plt.title('Full-spectrum Image')
        plt.subplot(1,2,2), plt.imshow(mymask), plt.title('Mask')
        if( [np.shape(Ttof)[0], np.shape(Ttof)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif(auto_mask):
        import skimage.filters
        mymask = reduction_tools.medianimage(Ttof)
        plt.figure()
        plt.subplot(1,3,1), plt.imshow(mymask), plt.title('Full-spectrum Image')
        mymask[mymask>mask_thresh[1]] = 0.0
        mymask[mymask<mask_thresh[0]] = 0.0
        mymask[mymask>0] = 1.0
        mymask[np.isinf(mymask)] = 0.0
        mymask[np.isnan(mymask)] = 0.0
        plt.subplot(1,3,2), plt.imshow(mymask), plt.title('Mask')
        mymask = skimage.filters.gaussian(mymask,sigma=2)
        mymask[mymask>0] = 1.0
        plt.subplot(1,3,3), plt.imshow(mymask), plt.title('Mask - gauss')
        plt.show(), plt.close()        
    else:
        mymask = np.ones([np.shape(Ttof)[0], np.shape(Ttof)[1]])
        
    if(debug_flag): #testing on a single pixel    
        T = Ttof[debug_idx[0],debug_idx[1],:]
        phasefitting_1D.phase_ratio_linearcomb(T,spectrum,phase1lac,phase2lac,phase_spectrum,lambda_range_norm,lambda_range_edges,est_phi=est_phi,method=method,bool_SG=bool_SG,SG_w=SG_w,SG_n=SG_n,bool_print=1)
        return

    phase_fraction = np.zeros(np.shape(mymask))
    for i in range(0, np.shape(mymask)[0]):
        print('---------------$$$$---------------')
        print('Processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                T = Ttof[i,j,:]
                try:
                    phi_fit = phasefitting_1D.phase_ratio_linearcomb(T,spectrum,phase1lac,phase2lac,phase_spectrum,lambda_range_norm,lambda_range_edges,est_phi=est_phi,method=method,bool_SG=bool_SG,SG_w=SG_w,SG_n=SG_n,bool_print=0)
                    phase_fraction[i,j] = phi_fit['phi']
                except:
                    print("Unexpected error at :", i, j)
                    phase_fraction[i,j] = -2.0

    if(bool_print):
        plt.imshow(phase_fraction), plt.title('Phase fraction (ph1 %)')
    if(bool_save):
        np.save('phase_fraction.npy', phase_fraction)
    return {'phase_fraction' : phase_fraction}

def phase_ratio_linearcomb_three_2D(Ttof,spectrum,phase1lac,phase2lac,phase3lac,phase_spectrum,lambda_range_norm,lambda_range_edges,filemask=0,auto_mask=True,mask_thresh=[0.05, 0.95],est_f1=0.3,est_f2=0.3,method='least_squares',bool_SG=False,SG_w=5,SG_n=1,bool_save=False,bool_print=False,debug_flag=False,debug_idx=[160,200]):
    """ Performs edge fitting with gaussian model to stack of TOF images (x,y,lambda)
    !!! IMPORTANT the data must be transmission or it won't work due to boundary conditions !!!
    
    INPUTS:
    Ttof = 3d matrix with the stack of TRANSMISSION (I/I0) TOF images (x,y,lambda) 
    spectrum_l = spectrum, length of this ndarray must correspond to size of Ttof(lambda)
    lambda_range = range corresponding to lambda where to perform the fitting
    filemask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_pos = estimated bragg edge position (in spectrum_l dimension)
    est_wid = estimated bragg edge width (in spectrum_l dimension)
    est_h = estimated bragg edge height (in spectrum_l dimension)
    bool_smooth = set to True to perform Savitzky-Golay filtering of the transmission derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_flag = set to True to test on single pixel (will actually also perform the whole image fitting if not stopped)
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    'edge_position' : edge position 
    'edge_height': edge height 
    'edge_width': edge width  
    'edge_slope': edge slope 
    'median_image': median Transmission image in the selected lambda range
    """ 

    if(filemask.any()):
        mymask = filemask
        plt.figure()
        plt.subplot(1,2,1), plt.imshow(np.median(Ttof,axis=2)), plt.title('Full-spectrum Image')
        plt.subplot(1,2,2), plt.imshow(mymask), plt.title('Mask')
        if( [np.shape(Ttof)[0], np.shape(Ttof)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif(auto_mask):
        import skimage.filters
        mymask = reduction_tools.medianimage(Ttof)
        plt.figure()
        plt.subplot(1,3,1), plt.imshow(mymask), plt.title('Full-spectrum Image')
        mymask[mymask>mask_thresh[1]] = 0.0
        mymask[mymask<mask_thresh[0]] = 0.0
        mymask[mymask>0] = 1.0
        mymask[np.isinf(mymask)] = 0.0
        mymask[np.isnan(mymask)] = 0.0
        plt.subplot(1,3,2), plt.imshow(mymask), plt.title('Mask')
        mymask = skimage.filters.gaussian(mymask,sigma=2)
        mymask[mymask>0] = 1.0
        plt.subplot(1,3,3), plt.imshow(mymask), plt.title('Mask - gauss')
        plt.show(), plt.close()        
    else:
        mymask = np.ones([np.shape(Ttof)[0], np.shape(Ttof)[1]])
        
    if(debug_flag): #testing on a single pixel    
        T = Ttof[debug_idx[0],debug_idx[1],:]
        phasefitting_1D.phase_ratio_linearcomb_three(T,spectrum,phase1lac,phase2lac,phase3lac,phase_spectrum,lambda_range_norm,lambda_range_edges,est_f1=est_f1,est_f2=est_f2,method=method,bool_SG=bool_SG,SG_w=SG_w,SG_n=SG_n,bool_print=1)
        return

    phase_fraction = np.zeros(np.shape(mymask))
    for i in range(0, np.shape(mymask)[0]):
        print('---------------$$$$---------------')
        print('Processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                T = Ttof[i,j,:]
                try:
                    phi_fit = phasefitting_1D.phase_ratio_linearcomb_three(T,spectrum,phase1lac,phase2lac,phase3lac,phase_spectrum,lambda_range_norm,lambda_range_edges,est_f1=est_f1,est_f2=est_f2,method=method,bool_SG=bool_SG,SG_w=SG_w,SG_n=SG_n,bool_print=0)
                    phase_fraction[i,j] = phi_fit['phi']
                except:
                    print("Unexpected error at :", i, j)
                    phase_fraction[i,j] = -2.0

    if(bool_print):
        plt.imshow(phase_fraction), plt.title('Phase fraction (ph1 %)')
    if(bool_save):
        np.save('phase_fraction.npy', phase_fraction)
    return {'phase_fraction' : phase_fraction}    