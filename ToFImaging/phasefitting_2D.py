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

from tqdm import tqdm
import time
   
def phase_ratio_linearcomb_2D(lac_tof,spectrum,phase1lac,phase2lac,spectrum_phase,lambda_range_norm,lambda_range_edges,calibration_matrix=np.ndarray([0]),mask=np.ndarray([0]),auto_mask=True,mask_thresh=[0.05, 0.95],est_phi=0.5,method='least_squares',bool_SG=False,SG_w=5,SG_n=1,bool_save=False,bool_print=False,debug_idx=[]):
    """ Performs phase ratio fitting on linear combination of two basis functions, works with linear attenuation coefficient (LAC) spectra
    INPUTS:
    lac_tof = 3d matrix with the stack of attenuation -log(I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
    phase1lac = lac of the phase 1
    phase2lac = lac of the phase 2
    spectrum_phase = spectrum corresponding to phase 1 and 2
    lambda_range_norm = lambda range where to normalize spectra
    lambda_range_edges = lambda range where to do the fitting
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y)
    mask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_phi = estimated phase fraction
    method = fitting method
    bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative
    SG_w = window size of S-G filter
    SG_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    dictionary with the following fit in the dimension of the mask
    'phase_ratio' : phase ratio
    """ 

    if(mask.any()):
        mymask = mask
        plt.figure()
        plt.subplot(1,2,1), plt.imshow(np.median(lac_tof,axis=2)), plt.title('Full-spectrum Image')
        plt.subplot(1,2,2), plt.imshow(mymask), plt.title('Mask')
        plt.tight_layout()
        plt.show()
        plt.close()        
        if( [np.shape(lac_tof)[0], np.shape(lac_tof)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif(auto_mask):
        import skimage.filters
        mymask = reduction_tools.median_image(lac_tof)
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
        plt.tight_layout()
        plt.show()
        plt.close()        
    else:
        mymask = np.ones([np.shape(lac_tof)[0], np.shape(lac_tof)[1]])
    
    if(calibration_matrix.any()):
        if ((np.shape(lac_tof)[0]!=np.shape(calibration_matrix)[0]) | (np.shape(lac_tof)[1]!=np.shape(calibration_matrix)[1])):
            print('!!!!!!! WARNING CALIBRATION MATRIX HAS NOT SAME SIZE OF IMAGE !!!!!!!!!!!!!!')
    
    if(any(debug_idx)): #testing on a single pixel    
        lac = lac_tof[debug_idx[0],debug_idx[1],:]

        if(calibration_matrix.any()):
            lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[debug_idx[0],debug_idx[1],1],calibration_matrix[debug_idx[0],debug_idx[1],0])
        else:
            lambd = spectrum

        phasefitting_1D.phase_ratio_linearcomb(lac,lambd,phase1lac,phase2lac,spectrum_phase,lambda_range_norm,lambda_range_edges,est_phi=est_phi,method=method,bool_SG=bool_SG,SG_w=SG_w,SG_n=SG_n,bool_print=1)
        return

    phase_ratio = np.zeros(np.shape(mymask))
    for i in tqdm(range(0, np.shape(mymask)[0])):
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                lac = lac_tof[i,j,:]

                if(calibration_matrix.any()):
                    lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[i,j,1],calibration_matrix[i,j,0])
                else:
                    lambd = spectrum

                try:
                    phi_fit = phasefitting_1D.phase_ratio_linearcomb(lac,lambd,phase1lac,phase2lac,spectrum_phase,lambda_range_norm,lambda_range_edges,est_phi=est_phi,method=method,bool_SG=bool_SG,SG_w=SG_w,SG_n=SG_n,bool_print=0)
                    phase_ratio[i,j] = phi_fit['phi']
                except:
                    print("Unexpected error at :", i, j)
                    phase_ratio[i,j] = -2.0

    if(bool_print):
        plt.figure()
        plt.imshow(phase_ratio, cmap='jet'), plt.title('Phase ratio (%)')
        plt.colorbar()
        plt.show(), plt.close()        
    if(bool_save):
        np.save('phase_ratio.npy', phase_ratio)
    return {'phase_ratio' : phase_ratio}

def phase_ratio_linearcomb_three_2D(lac_tof,spectrum,phase1lac,phase2lac,phase3lac,spectrum_phase,lambda_range_norm,lambda_range_edges,calibration_matrix=np.ndarray([0]),mask=np.ndarray([0]),auto_mask=True,mask_thresh=[0.05, 0.95],est_f1=0.333,est_f2=0.333,est_f3=0.334,method='least_squares',bool_SG=False,SG_w=5,SG_n=1,bool_save=False,bool_print=False,debug_idx=[]):
    """ Performs phase ratio fitting on linear combination of two basis functions, works with linear attenuation coefficient (LAC) spectra
    INPUTS:
    lac_tof = 3d matrix with the stack of attenuation -log(I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
    phase1lac = lac of the phase 1
    phase2lac = lac of the phase 2
    spectrum_phase = spectrum corresponding to phase 1 and 2
    lambda_range_norm = lambda range where to normalize spectra
    lambda_range_edges = lambda range where to do the fitting
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y)
    mask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_f1 = estimate phase 1 weight
    est_f2 = estimate phase 2 weight
    est_f3 = estimate phase 3 weight
    method = fitting method
    bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative
    SG_w = window size of S-G filter
    SG_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    dictionary with the following fit in the dimension of the mask
    'phase1_ratio' : phase1 weight
    'phase2_ratio' : phase2 weight
    'phase3_ratio' : phase3 weight
    """ 

    if(mask.any()):
        mymask = mask
        plt.figure()
        plt.subplot(1,2,1), plt.imshow(np.median(lac_tof,axis=2)), plt.title('Full-spectrum Image')
        plt.subplot(1,2,2), plt.imshow(mymask)
        plt.title('Mask')
        plt.tight_layout()
        plt.show()
        plt.close()
        if( [np.shape(lac_tof)[0], np.shape(lac_tof)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif(auto_mask):
        import skimage.filters
        mymask = reduction_tools.median_image(lac_tof)
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
        plt.tight_layout()
        plt.show()
        plt.close()
    else:
        mymask = np.ones([np.shape(lac_tof)[0], np.shape(lac_tof)[1]])
    
    if(calibration_matrix.any()):        
        if ((np.shape(lac_tof)[0]!=np.shape(calibration_matrix)[0]) | (np.shape(lac_tof)[1]!=np.shape(calibration_matrix)[1])):
            print('!!!!!!! WARNING CALIBRATION MATRIX HAS NOT SAME SIZE OF IMAGE !!!!!!!!!!!!!!')

    if(any(debug_idx)): #testing on a single pixel     
        lac = lac_tof[debug_idx[0],debug_idx[1],:]

        if(calibration_matrix.any()):     
            lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[debug_idx[0],debug_idx[1],1],calibration_matrix[debug_idx[0],debug_idx[1],0])
        else:
            lambd = spectrum

        phasefitting_1D.phase_ratio_linearcomb_three(lac,lambd,phase1lac,phase2lac,phase3lac,spectrum_phase,lambda_range_norm,lambda_range_edges,est_f1=est_f1,est_f2=est_f2,method=method,bool_SG=bool_SG,SG_w=SG_w,SG_n=SG_n,bool_print=1)
        return

    phase1_ratio = np.zeros(np.shape(mymask))
    phase2_ratio = np.zeros(np.shape(mymask))
    phase3_ratio = np.zeros(np.shape(mymask))
    for i in tqdm(range(0, np.shape(mymask)[0])):
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                lac = lac_tof[i,j,:]

                if(calibration_matrix.any()):     
                    lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[i,j,1],calibration_matrix[i,j,0])
                else:
                    lambd = spectrum

                try:
                    phi_fit = phasefitting_1D.phase_ratio_linearcomb_three(lac,lambd,phase1lac,phase2lac,phase3lac,spectrum_phase,lambda_range_norm,lambda_range_edges,est_f1=est_f1,est_f2=est_f2,method=method,bool_SG=bool_SG,SG_w=SG_w,SG_n=SG_n,bool_print=0)
                    phase1_ratio[i,j] = phi_fit['phi1']
                    phase2_ratio[i,j] = phi_fit['phi2']
                    phase3_ratio[i,j] = phi_fit['phi3']
                except:
                    print("Unexpected error at :", i, j)
                    phase1_ratio[i,j] = -2.0
                    phase2_ratio[i,j] = -2.0
                    phase3_ratio[i,j] = -2.0

    if(bool_print):
        plt.figure()
        plt.subplot(1,3,1), plt.imshow(phase1_ratio, cmap='jet'), plt.title('Phase 1 weight (%)'), plt.colorbar()
        plt.subplot(1,3,2), plt.imshow(phase2_ratio, cmap='jet'), plt.title('Phase 2 weight (%)'), plt.colorbar()
        plt.subplot(1,3,3), plt.imshow(phase3_ratio, cmap='jet'), plt.title('Phase 3 weight (%)'), plt.colorbar()
        plt.tight_layout()
        plt.show()
        plt.close()      
    if(bool_save):
        np.save('phase1_ratio.npy', phase1_ratio)
        np.save('phase2_ratio.npy', phase2_ratio)
        np.save('phase3_ratio.npy', phase3_ratio)
    return {'phase1_ratio' : phase1_ratio, 'phase2_ratio' : phase2_ratio, 'phase3_ratio' : phase3_ratio}    

def WavelengthSelectiveRatio2D(lac_tof,spectrum,l1,l2,l1_w=0,l2_w=0,calibration_matrix=np.ndarray([0]),mask=np.ndarray([0]),auto_mask=True,mask_thresh=[0.05, 0.95],bool_SG=False,SG_w=5,SG_n=1,bool_save=False,bool_print=False,debug_idx=[]):
    """ Performs phase ratio fitting on linear combination of two basis functions, works with linear attenuation coefficient (LAC) spectra
    INPUTS:
    lac_tof = 3d matrix with the stack of attenuation -log(I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
    l1 = lambda position of the 1st point
    l2 = lambda position of the 2nd point
    l1_w = lambda window size (bi-directional) of the 1st point
    l2_w = lambda window size (bi-directional) of the 2nd point
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y)
    mask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    dictionary with the following fit in the dimension of the mask
    'WSR' : Wavelength Selective Ratio
    """ 

    if(mask.any()):
        mymask = mask
        plt.figure()
        plt.subplot(1,2,1), plt.imshow(np.median(lac_tof,axis=2)), plt.title('Full-spectrum Image')
        plt.subplot(1,2,2), plt.imshow(mymask), plt.title('Mask')
        plt.tight_layout()
        plt.show()
        plt.close()        
        if( [np.shape(lac_tof)[0], np.shape(lac_tof)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif(auto_mask):
        import skimage.filters
        mymask = reduction_tools.median_image(lac_tof)
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
        plt.tight_layout()
        plt.show()
        plt.close()        
    else:
        mymask = np.ones([np.shape(lac_tof)[0], np.shape(lac_tof)[1]])
    
    if(calibration_matrix.any()):
        if ((np.shape(lac_tof)[0]!=np.shape(calibration_matrix)[0]) | (np.shape(lac_tof)[1]!=np.shape(calibration_matrix)[1])):
            print('!!!!!!! WARNING CALIBRATION MATRIX HAS NOT SAME SIZE OF IMAGE !!!!!!!!!!!!!!')
    
    if(any(debug_idx)): #testing on a single pixel    
        lac = lac_tof[debug_idx[0],debug_idx[1],:]

        if(calibration_matrix.any()):
            lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[debug_idx[0],debug_idx[1],1],calibration_matrix[debug_idx[0],debug_idx[1],0])
        else:
            lambd = spectrum

        phasefitting_1D.WavelengthSelectiveRatio(lac,lambd,l1,l2,l1_w,l2_w,bool_SG,SG_w,SG_n,bool_print=True)
        return

    phase_ratio = np.zeros(np.shape(mymask))
    for i in tqdm(range(0, np.shape(mymask)[0])):
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                lac = lac_tof[i,j,:]

                if(calibration_matrix.any()):
                    lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[i,j,1],calibration_matrix[i,j,0])
                else:
                    lambd = spectrum

                try:
                    phi_fit = phasefitting_1D.WavelengthSelectiveRatio(lac,lambd,l1,l2,l1_w,l2_w,bool_SG,SG_w,SG_n,bool_print=False)
                    phase_ratio[i,j] = phi_fit['WSR']
                except:
                    print("Unexpected error at :", i, j)
                    phase_ratio[i,j] = -2.0

    if(bool_print):
        plt.figure()
        plt.imshow(phase_ratio, cmap='jet'), plt.title('Phase ratio')
        plt.colorbar()
        plt.show(), plt.close()        
    if(bool_save):
        np.save('phase_ratio.npy', phase_ratio)
    return {'phase_ratio' : phase_ratio}    