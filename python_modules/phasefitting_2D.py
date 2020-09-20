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

import time
   
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
        mymask = reduction_tools.medianimage(LAC_tof)
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