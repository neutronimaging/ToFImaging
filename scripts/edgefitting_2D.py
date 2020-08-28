import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from skimage import io
import os, fnmatch
from os import listdir
from astropy.io import fits

import AdvancedBraggEdgeFitting
import TOF_routines
from TOF_routines import tof2l
from TOF_routines import find_nearest

import time
   
def image_edge_fitting_Tlambda(Ttof, spectrum_l, lambda_range, filemask=0, auto_mask = True, mask_thresh = [0.05, 0.95], est_pos=0, est_sigma=1, est_alpha=1, bool_save=True, bool_print=True, debug_flag=False, debug_idx=[160,200], bool_average=False, bool_linear=False):            
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
    if(filemask):
        mymask = io.imread(filemask)
        if( [np.shape(Ttof)[0], np.shape(Ttof)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif(auto_mask):
        import skimage.filters
        mymask = TOF_routines.medianimage(Ttof)
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
    myrange.append(find_nearest(spectrum_l, lambda_range[0])) # 3.7
    myrange.append(find_nearest(spectrum_l, lambda_range[1])) # 4.4
    if(est_pos):
        est_pos = find_nearest(spectrum_l[myrange[0]:myrange[1]], est_pos) # 4.05
    if(debug_flag): #testing on a single pixel    
        plt.imshow(Ttof.sum(axis=2))
        plt.show()
        plt.close()
        sp = np.zeros(len(spectrum_l))
        for i in range(0,len(spectrum_l)):
            sp[i] = np.median(Ttof[debug_idx[0],debug_idx[1],i]) # This is for the first Bragg edge fitting
        #run once the fitting to check if everything works
        AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(myspectrum=sp, myrange=myrange, myTOF=spectrum_l, est_pos=est_pos, est_sigma=est_sigma, est_alpha=est_alpha, bool_print=debug_flag, bool_average=bool_average, bool_linear=bool_linear)

    median_image = TOF_routines.medianimage(Ttof)
    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
        if(debug_flag):
            print('processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                print(i,j)
                mysignal = np.zeros(len(spectrum_l))
                for ind in range(0,len(spectrum_l)):
                    mysignal[ind] = np.median(Ttof[i,j,ind])
                try:
                    edge_fit = AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(myspectrum=mysignal, myrange=myrange, myTOF=spectrum_l, est_pos=est_pos, est_sigma=est_sigma, est_alpha=est_alpha, bool_print=False, bool_average=bool_average, bool_linear=bool_linear)
                    edge_position[i,j] = edge_fit['t0']
                    edge_height[i,j] = edge_fit['height']
                    if (len(edge_fit['pos_extrema'])==2):
                        edge_width[i,j] = spectrum_l[edge_fit['pos_extrema'][1]]-spectrum_l [edge_fit['pos_extrema'][0]]
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
        plt.imshow(edge_position)
        plt.figure()
        plt.imshow(edge_width)
        plt.figure()
        plt.imshow(edge_height)        
    if(bool_save):
        np.save('edge_position.npy', edge_position)
        np.save('edge_height.npy', edge_height)
        np.save('edge_width.npy', edge_width)
   
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width, 'median_image': median_image}        
    
def image_edge_fitting_Tlambda_gauss(Ttof, spectrum_l, lambda_range, filemask=0, auto_mask = True, mask_thresh = [0.05, 0.95], est_pos=0, est_wid=0, est_h=0, bool_smooth=False, smooth_w = 3, smooth_n = 0, bool_save=True, bool_print=True, debug_flag=False, debug_idx=[160,200]):        
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
    if(filemask):
        mymask = io.imread(filemask)
        if( [np.shape(Ttof)[0], np.shape(Ttof)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif(auto_mask):
        import skimage.filters
        mymask = TOF_routines.medianimage(Ttof)
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
    myrange.append(find_nearest(spectrum_l, lambda_range[0])) # 3.7
    myrange.append(find_nearest(spectrum_l, lambda_range[1])) # 4.4
    if(debug_flag): #testing on a single pixel    
        plt.imshow(Ttof.sum(axis=2))
        plt.show()
        plt.close()
        sp = Ttof[debug_idx[0],debug_idx[1],:]
        #run once the fitting to check if everything works
        AdvancedBraggEdgeFitting.GaussianBraggEdgeFitting(myspectrum=sp, myTOF=spectrum_l, myrange=myrange, est_pos = est_pos, est_wid=est_wid, est_h=est_h, bool_smooth=bool_smooth, smooth_w = smooth_w, smooth_n = smooth_n, bool_print=debug_flag)

    median_image = TOF_routines.medianimage(Ttof)
    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    edge_slope = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
        if(debug_flag):
            print('processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                #print(i,j)
                mysignal = Ttof[i,j,:]
                try:
                    edge_fit = AdvancedBraggEdgeFitting.GaussianBraggEdgeFitting(myspectrum=mysignal, myTOF=spectrum_l, myrange=myrange, est_pos = est_pos, est_wid=est_wid, est_h=est_h, bool_smooth=bool_smooth, smooth_w = smooth_w, smooth_n = smooth_n, bool_print=False)
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
        plt.figure(), plt.imshow(edge_position)
        plt.figure(), plt.imshow(edge_width)
        plt.figure(), plt.imshow(edge_height)        
        plt.figure(), plt.imshow(edge_slope)        
    if(bool_save):
        np.save('edge_position.npy', edge_position)
        np.save('edge_height.npy', edge_height)
        np.save('edge_width.npy', edge_width)
        np.save('edge_slope.npy', edge_slope)
   
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width, 'edge_slope': edge_slope, 'median_image': median_image}

def image_edge_fitting_T_gauss_calib(Ttof, spectrum, calibration_matrix, lambda_range, filemask=0, auto_mask = True, mask_thresh = [0.05, 0.95], est_pos=0, est_wid=0, est_h=0, bool_smooth=False, smooth_w = 3, smooth_n = 0, bool_save=True, bool_print=True, debug_flag=False, debug_idx=[160,200]):                
    """ Performs edge fitting with gaussian model to stack of TOF images (x,y,lambda)
    !!! IMPORTANT the data must be transmission or it won't work due to boundary conditions !!!
    
    INPUTS:
    Ttof = 3d matrix with the stack of TRANSMISSION (I/I0) TOF images (x,y,spectrum) 
    spectrum = spectrum, length of this ndarray must correspond to size of Ttof(lambda)
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum
    lambda_range = range corresponding to lambda where to perform the fitting
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
    #dictionary with the following fit in the dimension of the mask
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
        mymask = TOF_routines.medianimage(Ttof)
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
        plt.imshow(Ttof.sum(axis=2))
        plt.show()
        plt.close()
        lambd = TOF_routines.tof2l_calibration(spectrum,calibration_matrix[debug_idx[0],debug_idx[1],0],calibration_matrix[debug_idx[0],debug_idx[1],1])
        myrange = []
        myrange.append(find_nearest(lambd, lambda_range[0])) 
        myrange.append(find_nearest(lambd, lambda_range[1])) 
        
        T = Ttof[debug_idx[0],debug_idx[1],:] # This is for the first Bragg edge fitting
        #run once the fitting to check if everything works
        AdvancedBraggEdgeFitting.GaussianBraggEdgeFitting(myspectrum=T, myTOF=lambd, myrange=myrange, est_pos = est_pos, est_wid=est_wid, est_h=est_h, bool_smooth=bool_smooth, smooth_w = smooth_w, smooth_n = smooth_n, bool_print=debug_flag)

    median_image = TOF_routines.medianimage(Ttof)
    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    edge_slope = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
        if(debug_flag):
            print('processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                #print(i,j)
                T = Ttof[i,j,:]
                lambd = TOF_routines.tof2l_calibration(spectrum,calibration_matrix[debug_idx[0],debug_idx[1],0],calibration_matrix[debug_idx[0],debug_idx[1],1])
                myrange = []
                myrange.append(find_nearest(lambd, lambda_range[0])) 
                myrange.append(find_nearest(lambd, lambda_range[1])) 
                try:
                    edge_fit = AdvancedBraggEdgeFitting.GaussianBraggEdgeFitting(myspectrum=T, myTOF=spectrum_l, myrange=myrange, est_pos = est_pos, est_wid=est_wid, est_h=est_h, bool_smooth=bool_smooth, smooth_w = smooth_w, smooth_n = smooth_n, bool_print=False)
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
        plt.figure(), plt.imshow(edge_position)
        plt.figure(), plt.imshow(edge_width)
        plt.figure(), plt.imshow(edge_height)        
        plt.figure(), plt.imshow(edge_slope)        
    if(bool_save):
        np.save('edge_position.npy', edge_position)
        np.save('edge_height.npy', edge_height)
        np.save('edge_width.npy', edge_width)
        np.save('edge_slope.npy', edge_slope)
   
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width, 'edge_slope': edge_slope, 'median_image': median_image}

def image_phase_fraction(LAC_tof,spectrum_lambda,phase1lac,phase2lac,phase_spectrum,lambda_range_norm,lambda_range_edges,filemask=0, auto_mask = True, mask_thresh = [0.05, 0.95], debug_flag=0):        
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
    """ Performs FOBI reduction from open beam and sample spectrum    
    
    INPUTS:
    y = sample spectrum (I)
    y0 = open beam spectrum (I0)
    t = time-of-flight bins
    tmax = maximum time of flight (this parameter is dependent on the chopper frequency: tmax = 1/f  with f = chopper frequency)
    nrep = number of times the pattern is repeated (chopper)
    c = Wiener constant
    roll_value = shift the retrieved fobi spectra by this value 
    (this is necessary because the output from Wiener decorrelation places the peak of 
    the spectrum approximately in the middle of the array)
    SG_w = Savitzky Golay filter window
    SG_o = Savitzky Golay filter order

    OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    #'y_fobi': edge position
    #'x_fobi': edge height
    #'t_fobi': edge width
    """ 
    if(filemask):
        mymask = io.imread(filemask)
        if( [np.shape(LAC_tof)[0], np.shape(LAC_tof)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif(auto_mask):
        import skimage.filters
        mymask = TOF_routines.medianimage(LAC_tof)
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
        mymask = np.ones([np.shape(LAC_tof)[0], np.shape(LAC_tof)[1]])

    phase_fraction = np.zeros(np.shape(mymask))
    for i in range(0, np.shape(mymask)[0]):
            if(debug_flag):
                print('processing row n. ', i, 'of', np.shape(mymask)[0])
            for j in range(0, np.shape(mymask)[1]):
                if (mymask[i,j]):
                    lac = LAC_tof[i,j,:]
                    try:
                        phi_fit = AdvancedBraggEdgeFitting.phase_fraction_fitting(lac,spectrum_lambda,phase1lac,phase2lac,phase_spectrum,lambda_range_norm,lambda_range_edges,bool_plot=0)
                        phase_fraction[i,j] = phi_fit['phi']
                    except:
                        print("Unexpected error at :", i, j)
                        phase_fraction[i,j] = -2.0
    return {'phase_fraction' : phase_fraction}