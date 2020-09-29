import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from skimage import io
import os, fnmatch
from os import listdir
from astropy.io import fits

import edgefitting_1D
import reduction_tools
from reduction_tools import tof2l
from reduction_tools import find_nearest

from tqdm import tqdm
import time
   
def AdvancedBraggEdgeFitting_2D(Ttof,spectrum,spectrum_range,calibration_matrix=np.ndarray([0]),mask=np.ndarray([0]),auto_mask=True,mask_thresh=[0.05, 0.95],est_pos=0,est_sigma=1,est_alpha=1,bool_smooth=True,smooth_w=5,smooth_n=1,bool_linear=False,bool_save=False,bool_print=False,debug_idx=[]):            
    """ Performs edge fitting with gaussian model to stack of TOF images (x,y,lambda)
    INPUTS:
    Ttof = 3d matrix with the stack of TRANSMISSION (I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of Ttof(lambda)
    spectrum_range = range corresponding to lambda where to perform the fitting
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y)
    mask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_pos = estimated bragg edge position (in spectrum dimension)
    est_sigma = expected Gaussian broadening
    est_alpha = expected decay constant (moderator property)
    bool_smooth = set to True to perform Savitzky-Golay filtering of the transmission derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    Dictionary with the following fit in the dimension of the mask
    'edge_position' : edge position 
    'edge_height': edge height 
    'edge_width': edge width  
    'edge_slope': edge slope 
    'median_image': median Transmission image in the selected lambda range
    """ 

    if(mask.any()):
        mymask = mask
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

    if(calibration_matrix.any()):
        if ((np.shape(Ttof)[0]!=np.shape(calibration_matrix)[0]) | (np.shape(Ttof)[1]!=np.shape(calibration_matrix)[1])):
            print('!!!!!!! WARNING CALIBRATION MATRIX HAS NOT SAME SIZE OF IMAGE !!!!!!!!!!!!!!')

    if(any(debug_idx)): #testing on a single pixel  
        signal = Ttof[debug_idx[0],debug_idx[1],:]

        if(calibration_matrix.any()):
            lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[debug_idx[0],debug_idx[1],1],calibration_matrix[debug_idx[0],debug_idx[1],0])
        else:
            lambd = spectrum

        edgefitting_1D.AdvancedBraggEdgeFitting(signal=signal,spectrum=lambd,spectrum_range=spectrum_range,est_pos=est_pos,est_sigma=est_sigma,est_alpha=est_alpha,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_linear=bool_linear,bool_print=True)
        return

    median_image = reduction_tools.medianimage(Ttof)
    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    
    start_time = time.time()
    for i in tqdm(range(0, np.shape(mymask)[0])): #loop for all pixel position, where the mask is equal to one
        # print('---------------$$$$---------------')
        # print('Processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                signal = Ttof[i,j,:]

                if(calibration_matrix.any()):
                    lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[i,j,1],calibration_matrix[i,j,0])
                else:
                    lambd = spectrum

                try:
                    edge_fit = edgefitting_1D.AdvancedBraggEdgeFitting(signal=signal,spectrum=lambd,spectrum_range=spectrum_range,est_pos=est_pos,est_sigma=est_sigma,est_alpha=est_alpha,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_linear=bool_linear,bool_print=False)
                    edge_position[i,j] = edge_fit['t0']
                    edge_height[i,j] = edge_fit['height']
                    if (len(edge_fit['pos_extrema'])==2):
                        edge_width[i,j] = spectrum[edge_fit['pos_extrema'][1]]-spectrum [edge_fit['pos_extrema'][0]]
                    else:
                        edge_width[i,j] = -2.0
                except:
                    print("Unexpected error at :", i, j)
                    edge_position[i,j] = -2.0
                    edge_width[i,j] = -2.0
                    edge_height[i,j] = -2.0
    print("--- %s seconds ---" % (time.time() - start_time))

    if(bool_print):
        plt.figure()
        plt.subplot(1,3,1), plt.imshow(edge_position), plt.title('Edge position')
        plt.subplot(1,3,2), plt.imshow(edge_width), plt.title('Edge width')
        plt.subplot(1,3,3), plt.imshow(edge_height), plt.title('Edge height')
        plt.show(), plt.close()        
    if(bool_save):
        np.save('edge_position.npy', edge_position)
        np.save('edge_height.npy', edge_height)
        np.save('edge_width.npy', edge_width)
   
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width, 'median_image': median_image}    

def GaussianBraggEdgeFitting_2D(Ttof,spectrum,spectrum_range,calibration_matrix=np.ndarray([0]),mask=np.ndarray([0]),auto_mask=True,mask_thresh=[0.05, 0.95],est_pos=0,est_wid=0,est_h=0,pos_BC=0,wid_BC=0,h_BC=0,bool_log=False,bool_smooth=True,smooth_w=5,smooth_n=1,bool_save=False,bool_print=False,debug_idx=[]):        
    """ Performs edge fitting with gaussian model to stack of TOF images (x,y,lambda)
    INPUTS:
    Ttof = 3d matrix with the stack of TRANSMISSION (I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of Ttof(lambda)
    spectrum_range = range corresponding to lambda where to perform the fitting
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y)
    mask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_pos = estimated bragg edge position (in spectrum dimension)
    est_wid = estimated bragg edge width (in spectrum dimension)
    est_h = estimated bragg edge height (in spectrum dimension)
    pos_BC = boundary conditions for the bragg edge position fit
    wid_BC = boundary conditions for the bragg edge width fit
    h_BC = boundary conditions for the bragg edge height fit
    bool_log = set to True to convert T to mu
    bool_smooth = set to True to perform Savitzky-Golay filtering of the transmission derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    Dictionary with the following fit in the dimension of the mask
    'edge_position' : edge position 
    'edge_height': edge height 
    'edge_width': edge width  
    'edge_slope': edge slope 
    'median_image': median Transmission image in the selected lambda range
    """ 

    if(mask.any()):
        mymask = mask
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

    if(bool_log):
        Ttof = -np.log(Ttof)

    if(calibration_matrix.any()):
        if ((np.shape(Ttof)[0]!=np.shape(calibration_matrix)[0]) | (np.shape(Ttof)[1]!=np.shape(calibration_matrix)[1])):
            print('!!!!!!! WARNING CALIBRATION MATRIX HAS NOT SAME SIZE OF IMAGE !!!!!!!!!!!!!!')

    if(any(debug_idx)): #testing on a single pixel    
        signal = Ttof[debug_idx[0],debug_idx[1],:]

        if(calibration_matrix.any()):
            lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[debug_idx[0],debug_idx[1],1],calibration_matrix[debug_idx[0],debug_idx[1],0])
        else:
            lambd = spectrum

        edgefitting_1D.GaussianBraggEdgeFitting(signal=signal,spectrum=lambd,spectrum_range=spectrum_range,est_pos=est_pos,est_wid=est_wid,est_h=est_h,pos_BC=pos_BC,wid_BC=wid_BC,h_BC=h_BC,bool_log=bool_log,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_print=True)
        return

    median_image = reduction_tools.medianimage(Ttof)
    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    edge_slope = np.zeros(np.shape(mymask))
    
    start_time = time.time()
    for i in tqdm(range(0, np.shape(mymask)[0])): #loop for all pixel position, where the mask is equal to one
        # print('---------------$$$$---------------')
        # print('Processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                signal = Ttof[i,j,:]

                if(calibration_matrix.any()):
                    lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[i,j,1],calibration_matrix[i,j,0])
                else:
                    lambd = spectrum

                try:
                    edge_fit = edgefitting_1D.GaussianBraggEdgeFitting(signal=signal,spectrum=lambd,spectrum_range=spectrum_range,est_pos=est_pos,est_wid=est_wid,est_h=est_h,pos_BC=pos_BC,wid_BC=wid_BC,h_BC=h_BC,bool_log=bool_log,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_print=False)
                    edge_position[i,j] = edge_fit['t0']
                    edge_height[i,j] = edge_fit['edge_height']
                    edge_width[i,j] = edge_fit['edge_width']
                    edge_slope[i,j] = edge_fit['edge_slope']
                except:
                    print("Unexpected error at :", i, j)
                    edge_position[i,j] = -2.0
                    edge_width[i,j] = -2.0
                    edge_height[i,j] = -2.0
                    edge_slope[i,j] = -2.0
    print("--- %s seconds ---" % (time.time() - start_time))

    if(bool_print):
        plt.figure()
        plt.subplot(1,3,1), plt.imshow(edge_position), plt.title('Edge position')
        plt.subplot(1,3,2), plt.imshow(edge_width), plt.title('Edge width')
        plt.subplot(1,3,3), plt.imshow(edge_height), plt.title('Edge height')
        plt.show(), plt.close()       
    if(bool_save):
        np.save('edge_position.npy', edge_position)
        np.save('edge_height.npy', edge_height)
        np.save('edge_width.npy', edge_width)
        np.save('edge_slope.npy', edge_slope)
   
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width, 'edge_slope': edge_slope, 'median_image': median_image}
