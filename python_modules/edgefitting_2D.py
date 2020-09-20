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

import time
   
def AdvancedBraggEdgeFitting_2D(Ttof,spectrum,spectrum_range,filemask=0,auto_mask=True,mask_thresh=[0.05, 0.95],est_pos=0,est_sigma=1,est_alpha=1,bool_smooth=True,smooth_w=5,smooth_n=1,bool_linear=False,bool_save=False,bool_print=False,debug_flag=False,debug_idx=[160,200]):            
    """ Performs edge fitting with gaussian model to stack of TOF images (x,y,lambda)
    
    INPUTS:
    Ttof = 3d matrix with the stack of TRANSMISSION (I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of Ttof(lambda)
    spectrum_range = range corresponding to lambda where to perform the fitting
    filemask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_pos = estimated bragg edge position (in spectrum dimension)
    est_wid = estimated bragg edge width (in spectrum dimension)
    est_h = estimated bragg edge height (in spectrum dimension)
    bool_smooth = set to True to perform Savitzky-Golay filtering of the transmission derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_flag = set to True to test on single pixel (will actually also perform the whole image fitting if not stopped)
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    Dictionary with the following fit in the dimension of the mask
    'edge_position' : edge position 
    'edge_height': edge height 
    'edge_width': edge width  
    'edge_slope': edge slope 
    'median_image': median Transmission image in the selected lambda range
    """ 
    if(filemask):
        mymask = io.imread(filemask)
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

    myrange = []
    myrange.append(find_nearest(spectrum, spectrum_range[0])) # 3.7
    myrange.append(find_nearest(spectrum, spectrum_range[1])) # 4.4
    if(est_pos):
        est_pos = find_nearest(spectrum[myrange[0]:myrange[1]], est_pos) # 4.05
    if(debug_flag): #testing on a single pixel  
        signal = Ttof[debug_idx[0],debug_idx[1],:]
        try:
            edgefitting_1D.AdvancedBraggEdgeFitting(signal=signal,spectrum=spectrum,spectrum_range=myrange,est_pos=est_pos,est_sigma=est_sigma,est_alpha=est_alpha,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_linear=bool_linear,bool_print=True)
        except:
            print('Fitting failed')
        return

    median_image = reduction_tools.medianimage(Ttof)
    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
        print('---------------$$$$---------------')
        print('Processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                signal = Ttof[i,j,:]
                try:
                    edge_fit = edgefitting_1D.AdvancedBraggEdgeFitting(signal=signal,spectrum=spectrum,spectrum_range=myrange,est_pos=est_pos,est_sigma=est_sigma,est_alpha=est_alpha,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_linear=bool_linear,bool_print=False)
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
    
def GaussianBraggEdgeFitting_2D(Ttof,spectrum,spectrum_range,filemask=0,auto_mask=True,mask_thresh=[0.05, 0.95],est_pos=0,est_wid=0,est_h=0,bool_log=False,bool_smooth=True,smooth_w=5,smooth_n=1,bool_save=False,bool_print=False,debug_flag=False,debug_idx=[160,200]):        
    """ Performs edge fitting with gaussian model to stack of TOF images (x,y,lambda)
    
    INPUTS:
    Ttof = 3d matrix with the stack of TRANSMISSION (I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of Ttof(lambda)
    spectrum_range = range corresponding to lambda where to perform the fitting
    filemask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_pos = estimated bragg edge position (in spectrum dimension)
    est_wid = estimated bragg edge width (in spectrum dimension)
    est_h = estimated bragg edge height (in spectrum dimension)
    bool_smooth = set to True to perform Savitzky-Golay filtering of the transmission derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    bool_log = set to True to convert T to mu
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_flag = set to True to test on single pixel (will actually also perform the whole image fitting if not stopped)
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    Dictionary with the following fit in the dimension of the mask
    'edge_position' : edge position 
    'edge_height': edge height 
    'edge_width': edge width  
    'edge_slope': edge slope 
    'median_image': median Transmission image in the selected lambda range
    """ 
    if(filemask):
        mymask = io.imread(filemask)
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

    myrange = []
    myrange.append(find_nearest(spectrum, spectrum_range[0])) # 3.7
    myrange.append(find_nearest(spectrum, spectrum_range[1])) # 4.4
    if(debug_flag): #testing on a single pixel    
        sp = Ttof[debug_idx[0],debug_idx[1],:]
        try:
            edgefitting_1D.GaussianBraggEdgeFitting(signal=sp,spectrum=spectrum,spectrum_range=myrange,est_pos=est_pos,est_wid=est_wid,est_h=est_h,bool_log=bool_log,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_print=debug_flag)
        except:
            print('Fitting failed')
        return

    median_image = reduction_tools.medianimage(Ttof)
    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    edge_slope = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
        print('---------------$$$$---------------')
        print('Processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                sp = Ttof[i,j,:]
                try:
                    edge_fit = edgefitting_1D.GaussianBraggEdgeFitting(signal=sp,spectrum=spectrum,spectrum_range=myrange,est_pos=est_pos,est_wid=est_wid,est_h=est_h,bool_log=bool_log,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_print=False)
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

def GaussianBraggEdgeFitting_2D_Calib_matrix(Ttof,spectrum,calibration_matrix,spectrum_range,filemask=0,auto_mask=True,mask_thresh=[0.05, 0.95],est_pos=0,est_wid=0,est_h=0,bool_log=False,bool_smooth=True,smooth_w=5,smooth_n=1,bool_save=False,bool_print=False,debug_flag=False,debug_idx=[160,200]):        
    """ Performs edge fitting with gaussian model to stack of TOF images (x,y,lambda)
    
    INPUTS:
    Ttof = 3d matrix with the stack of TRANSMISSION (I/I0) TOF images (x,y,spectrum) 
    spectrum = spectrum, length of this ndarray must correspond to size of Ttof(lambda)
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y)
    spectrum_range = range corresponding to lambda where to perform the fitting
    filemask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_pos = estimated bragg edge position (in spectrum dimension)
    est_wid = estimated bragg edge width (in spectrum dimension)
    est_h = estimated bragg edge height (in spectrum dimension)
    bool_smooth = set to True to perform Savitzky-Golay filtering of the transmission derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_flag = set to True to test on single pixel (will actually also perform the whole image fitting if not stopped)
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    Dictionary with the following fit in the dimension of the mask
    'edge_position' : edge position 
    'edge_height': edge height 
    'edge_width': edge width  
    'edge_slope': edge slope 
    'median_image': median Transmission image in the selected lambda range
    """ 

    if(filemask):
        mymask = io.imread(filemask)
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

    if (np.shape(Ttof)[0]!=np.shape(calibration_matrix)[0] | np.shape(Ttof)[1]!=np.shape(calibration_matrix)[1]):
        print('!!!!!!! WARNING CALIBRATION MATRIX HAS NOT SAME SIZE OF IMAGE !!!!!!!!!!!!!!')
    if(debug_flag): #testing on a single pixel    
        lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[debug_idx[0],debug_idx[1],1],calibration_matrix[debug_idx[0],debug_idx[1],0])
        myrange = []
        myrange.append(find_nearest(lambd, spectrum_range[0])) 
        myrange.append(find_nearest(lambd, spectrum_range[1])) 
        
        signal = Ttof[debug_idx[0],debug_idx[1],:]
        try:
            edgefitting_1D.GaussianBraggEdgeFitting(signal=sp,spectrum=spectrum,spectrum_range=myrange,est_pos=est_pos,est_wid=est_wid,est_h=est_h,bool_log=bool_log,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_print=debug_flag)
        except:
            print('Fitting failed')
        return


    median_image = reduction_tools.medianimage(Ttof)
    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    edge_slope = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
        print('---------------$$$$---------------')
        print('Processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                signal = Ttof[i,j,:]
                lambd = reduction_tools.tof2l_t0k(spectrum,calibration_matrix[i,j,1],calibration_matrix[i,j,0])
                myrange = []
                myrange.append(find_nearest(lambd, spectrum_range[0])) 
                myrange.append(find_nearest(lambd, spectrum_range[1])) 
                try:
                    edge_fit = edgefitting_1D.GaussianBraggEdgeFitting(signal=signal,spectrum=spectrum,spectrum_range=myrange,est_pos=est_pos,est_wid=est_wid,est_h=est_h,bool_log=bool_log,bool_smooth=bool_smooth,smooth_w=smooth_w,smooth_n=smooth_n,bool_print=False)
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
