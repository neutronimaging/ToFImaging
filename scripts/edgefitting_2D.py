
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import AdvancedBraggEdgeFitting
# from PIL import Image

import tifffile
from tifffile import TiffFile

from skimage import io

import os, fnmatch
from os import listdir

from TOF_routines import tof2l
from TOF_routines import find_nearest

import time


# #path to sample images
# pathdata = '/media/carminati_c/Data2/LSP_Manuel/sample_binned_spotCleaned/'

# #path to open beam images
# pathob = '/media/carminati_c/Data2/LSP_Manuel/OB_binned_spotCleaned/'

# pathspectrum = '/media/ws_niag/10_people/Morgano/RADEN_data_analysis/TEST6_000_Spectra.txt'

# pathmask = '/pathtomask/'
# filemask = pathmask+'mask.tif'

def image_edge_fitting(pathdata, pathob, filemask, pathspectrum, cal_parameters, binning, lambda_range, initial_guess, bool_save):
    """ Run edge_fitting on provided 3D image, the third dimension is supposed to be lambda(A) or TOF 
    inputs:
    pathdata: path to the sample data to be processed
    pathob: path to the open beam images
    filemask: path to the image to be used as mask for fitting
    pathspectrum: path the spectrum ile
    cal_paramters: dictionary containing the calibration parameters, t0 (delay) and L (flight path)
    binning: binning factor to be applied to the spectrum file, images are supposed to be already binned
    lambda_range: 2 values array with first and last lambda value to be considered, defines the area of the spectrum to be analysed
    initial_guess: dictionary containing the initial guess for the edge parameters: to (edge position), sigma (gaussian broadening), tau (moderator decay)
    bool_save: boolean value that triggers the saving of the results
    
    returns: 2D map for edge position, edge height and edge width
    """
    
    files_sample = (sorted(fnmatch.filter(listdir(pathdata),'*.tif')))
    files_ob = (sorted(fnmatch.filter(listdir(pathob),'*.tif')))
    mymask = io.imread(filemask)

    #load the spectrum
    spectrum = np.loadtxt(pathspectrum, usecols=0)


    spectrum_binned = spectrum[0::binning]

    #here che calibration parameters
#     t0=-0.0002618673892937752
#     L=18.961609065251505
    
    t0 = cal_parameters.get('t0')
    L = cal_parameters.get('L')
    mylambda_bin = tof2l(spectrum_binned,t0,L)

    ob = io.imread(pathob+'/'+files_ob[0])

    ob_image = np.zeros([np.shape(ob)[0], np.shape(ob)[1], len(files_ob)])
    sample_image = np.zeros([np.shape(ob)[0], np.shape(ob)[1], len(files_ob)])


    #read the images, in this example they are already selected in a neighborhood of the studied edge 
    for i in range(0,len(files_ob)):
        ob_image[:,:,i] = io.imread(pathob+'/'+files_ob[i])
        sample_image[:,:,i] = io.imread(pathdata+'/'+files_sample[i])

    #compute transmission image
    trans_image = np.zeros([np.shape(ob)[0], np.shape(ob)[1],len(files_ob)])
    for i in range(0,len(files_ob)):
        dose = np.median(ob_image[10:20,10:20,i])/np.median(sample_image[10:20,10:20,i]) # update the doseROI
        trans_image[:,:,i] = sample_image[:,:,i]/ob_image[:,:,i]*dose

    sp=np.zeros(len(files_ob))
    for i in range(0,len(files_ob)):
        sp[i] = np.median(trans_image[44,204,i]) # This is for the first Bragg edge fitting


    #here define the lambda range to be used for single edge fitting
    myrange =[]
    myrange.append(find_nearest(mylambda_bin, lambda_range[0])) # 3.7
    myrange.append(find_nearest(mylambda_bin, lambda_range[1])) # 4.4
  


    est_sigma = initial_guess.get('sigma')#0.08
    est_alpha = initial_guess.get('alpha')#0.08
    est_pos = find_nearest(mylambda_bin[myrange[0]:myrange[1]], initial_guess.get('t0')) # 4.05
    small_range = np.array([0, myrange[1]-myrange[0]-1])
    small_lambda = mylambda_bin[myrange[0]:myrange[1]]

    #run once the fitting to check if everything works
    AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(sp[myrange[0]:myrange[1]], small_range, small_lambda, est_pos, est_sigma, est_alpha, False, False, False, True)

    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_heigth = np.zeros(np.shape(mymask))

    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):

#         print('processing row n. ', i, 'of', np.shape(mymask)[0])

        for j in range(0, np.shape(mymask)[1]):

    #         if (little_mask[i,j]):
            if (mymask[i,j]):
    #             print(i,j,' ciao')
                # extract the signal
                mysignal = np.zeros(myrange[1]-myrange[0])


                for ind in range(myrange[0],myrange[1]):
                    mysignal[ind-myrange[0]] = np.median(trans_image[i,j,ind])

                try:
                    edge_fit = AdvancedBraggEdgeFitting_v2.AdvancedBraggEdgeFitting(mysignal, small_range, small_lambda, est_pos, est_sigma, est_alpha, False, False, False, True)

                    edge_position[i,j] = edge_fit['t0']
                    edge_heigth[i,j] = edge_fit['height']
                    if (len(edge_fit['pos_extrema'])==2):
                        edge_width[i,j] = small_lambda[edge_fit['pos_extrema'][1]]-small_lambda[edge_fit['pos_extrema'][0]]
                    else:
                        edge_width[i,j] = -2.0
                except:
                    print("Unexpected error at :", i, j)
                    edge_position[i,j] = -2.0
                    edge_width[i,j] = -2.0
                    edge_heigth[i,j] = -2.0


    # plt.figure()
    # plt.imshow(edge_position)
    # save the results
    if (bool_save):
        np.save('edge_position.npy', edge_position)
        np.save('edge_height.npy', edge_heigth)
        np.save('edge_width.npy', edge_width)


    print("--- %s seconds ---" % (time.time() - start_time))
    
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width}

