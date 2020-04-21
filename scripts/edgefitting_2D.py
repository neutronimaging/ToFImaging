import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import tifffile
from tifffile import TiffFile
from skimage import io
import os, fnmatch
from os import listdir
from astropy.io import fits

import AdvancedBraggEdgeFitting
# from PIL import Image
from TOF_routines import tof2l
from TOF_routines import find_nearest

import time

# #path to sample images
# pathdata = '/media/carminati_c/Data2/LSP_Manuel/sample_binned_spotCleaned/'

# #path to open beam images
# pathob = '/media/carminati_c/Data2/LSP_Manuel/OB_binned_spotCleaned/'

# pathspectrum = '/media/ws_niag/10_people/Morgano/RADEN_data_analysis/TEST6_000_Spectra.txt'

# pathmask = '/pathmask/'
# filemask = pathmask+'mask.tif'

def image_edge_fitting(pathdata, pathob, pathspectrum, lambda_range, filemask=0, t0_cal=-0.0002618673892937752, L_cal=18.961609065251505, binning=1, est_pos=0, est_sigma=1, est_alpha=1, bool_dose=0, bool_save=0, bool_print=1):
    ##INPUTS:
    #pathdata: path to sample images
    #pathob: path to open beam images 
    #filemask: path to the .tif mask of the region to fit. if not specified runs over the whole area
    #pathspectrum: path to the .txt spectrum
    #t0_cal: t0 for conversion from tof to lambda
    #L_cal: L for conversion from tof to lambda
    #binning: binning periodicity
    #lambda_range: lambda range for single edge fitting (in Angstrom): 2 values array
    #est_pos: expected bragg edge position (in A), if unset is found automatically
    #est_sigma: expected Gaussian broadening
    #est_alpha: expected decay constant (moderator property))
    #bool_dose: flag for dose correction
    #bool_save: flag for data saving
    #bool_print: flag for image printing
    
    ##OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    #'edge_position' : edge position
    #'edge_height': edge height
    #'edge_width': edge width
    
    debug_flag = 1
    
    #files_sample = (sorted(fnmatch.filter(listdir(pathdata),'*.tif')))
    #files_ob = (sorted(fnmatch.filter(listdir(pathob),'*.tif')))
    files_sample = (sorted(fnmatch.filter(listdir(pathdata),'*.fits')))
    files_ob = (sorted(fnmatch.filter(listdir(pathob),'*.fits')))
    
    #load the spectrum
    spectrum = np.loadtxt(pathspectrum, usecols=0)
    spectrum_binned = spectrum[0::binning]
    
    mylambda_bin = tof2l(spectrum_binned,t0_cal,L_cal)

    #ob = io.imread(pathob+'/'+files_ob[0])
    ob = fits.getdata(pathob+'/'+files_ob[0])
    
    ob_image = np.zeros([np.shape(ob)[0], np.shape(ob)[1], len(files_ob)])
    sample_image = np.zeros([np.shape(ob)[0], np.shape(ob)[1], len(files_ob)])
    
    if(filemask):
        mymask = io.imread(filemask)
    else:
        mymask = np.ones([np.shape(ob)[0], np.shape(ob)[1], len(files_ob)])

    #read the images, in this example they are already selected in a neighborhood of the studied edge 
    for i in range(0,len(files_ob)):
        #ob_image[:,:,i] = io.imread(pathob+'/'+files_ob[i])
        #sample_image[:,:,i] = io.imread(pathdata+'/'+files_sample[i])
        ob_image[:,:,i] = fits.getdata(pathob+'/'+files_ob[i])
        sample_image[:,:,i] = fits.getdata(pathdata+'/'+files_sample[i])

    #compute transmission image
    trans_image = np.zeros([np.shape(ob)[0], np.shape(ob)[1],len(files_ob)])
    dose = 1
    for i in range(0,len(files_ob)):
        if(bool_dose):
            dose = np.median(ob_image[10:20,10:20,i])/np.median(sample_image[10:20,10:20,i]) # update the doseROI
        trans_image[:,:,i] = sample_image[:,:,i]/ob_image[:,:,i]*dose

    #here define the lambda range to be used for single edge fitting
    myrange = []
    myrange.append(find_nearest(mylambda_bin, lambda_range[0])) # 3.7
    myrange.append(find_nearest(mylambda_bin, lambda_range[1])) # 4.4
    if(est_pos):
        est_pos = find_nearest(mylambda_bin[myrange[0]:myrange[1]], est_pos) # 4.05
    
    if(debug_flag): #testing on a single pixel    
        plt.imshow(trans_image.sum(axis=2))
        plt.show()
        plt.close()
        small_range = np.array([0, myrange[1]-myrange[0]-1])
        small_lambda = mylambda_bin[myrange[0]:myrange[1]]
        sp = np.zeros(len(files_ob))
        for i in range(0,len(files_ob)):
            sp[i] = np.median(trans_image[50,50,i]) # This is for the first Bragg edge fitting
        #run once the fitting to check if everything works
        AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(myspectrum=sp[myrange[0]:myrange[1]], myrange=small_range, myTOF=small_lambda, est_pos=est_pos, est_sigma=est_sigma, est_alpha=est_alpha, bool_print=1, bool_average=0, bool_linear=0)

    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_heigth = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
#         print('processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
    #             print(i,j,' ciao')
                # extract the signal
                #mysignal = np.zeros(myrange[1]-myrange[0])
                mysignal = np.zeros(np.int(np.shape(ob)[2]))
                #for ind in range(myrange[0],myrange[1]):
                    #mysignal[ind-myrange[0]] = np.median(trans_image[i,j,ind])
                for ind in range(0,np.int(np.shape(mysignal)[0])):
                    mysignal[ind] = np.median(trans_image[i,j,ind])
                try:
                    edge_fit = AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(myspectrum=mysignal, myrange=myrange, myTOF=mylambda_bin, est_pos=est_pos, est_sigma=est_sigma, est_alpha=est_alpha, bool_print=0, bool_average=0, bool_linear=0)
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

    if(bool_print):
        plt.figure()
        plt.imshow(edge_position)
        plt.figure()
        plt.imshow(edge_width)
        plt.figure()
        plt.imshow(edge_heigth)        
    if(bool_save):
        np.save('edge_position.npy', edge_position)
        np.save('edge_height.npy', edge_heigth)
        np.save('edge_width.npy', edge_width)

    print("--- %s seconds ---" % (time.time() - start_time))
   
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width}

def image_edge_fitting_T(pathdata, pathspectrum, lambda_range, filemask=0, t0_cal=-0.0002618673892937752, L_cal=18.961609065251505, est_pos=0, est_sigma=1, est_alpha=1, bool_save=0, bool_print=1):
    ##INPUTS:
    #pathdata: path to sample images
    #pathob: path to open beam images 
    #filemask: path to the .tif mask of the region to fit. if not specified runs over the whole area
    #pathspectrum: path to the .txt spectrum
    #t0_cal: t0 for conversion from tof to lambda
    #L_cal: L for conversion from tof to lambda
    #binning: binning periodicity
    #lambda_range: lambda range for single edge fitting (in Angstrom): 2 values array
    #est_pos: expected bragg edge position (in A), if unset is found automatically
    #est_sigma: expected Gaussian broadening
    #est_alpha: expected decay constant (moderator property))
    #bool_dose: flag for dose correction
    #bool_save: flag for data saving
    #bool_print: flag for image printing
    
    ##OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    #'edge_position' : edge position
    #'edge_height': edge height
    #'edge_width': edge width
    
    debug_flag = 1
    
    files_sample = (sorted(fnmatch.filter(listdir(pathdata),'*.tif')))
    
    #load the spectrum
    spectrum = np.loadtxt(pathspectrum, usecols=0)
    
    mylambda_bin = tof2l(spectrum,t0_cal,L_cal)

    ob = io.imread(pathdata+'/'+files_sample[0])
    
    sample_image = np.zeros([np.shape(ob)[0], np.shape(ob)[1], len(files_sample)])
    
    if(filemask):
        mymask = io.imread(filemask)
    else:
        mymask = np.ones([np.shape(ob)[0], np.shape(ob)[1]])

    #read the images, in this example they are already selected in a neighborhood of the studied edge 
    for i in range(0,len(files_sample)):
        sample_image[:,:,i] = io.imread(pathdata+'/'+files_sample[i])/256

    #here define the lambda range to be used for single edge fitting
    myrange = []
    myrange.append(find_nearest(mylambda_bin, lambda_range[0])) # 3.7
    myrange.append(find_nearest(mylambda_bin, lambda_range[1])) # 4.4
    if(est_pos):
        est_pos = find_nearest(mylambda_bin[myrange[0]:myrange[1]], est_pos) # 4.05
    
    if(debug_flag): #testing on a single pixel    
        plt.imshow(sample_image.sum(axis=2))
        plt.show()
        plt.close()
        sp = np.zeros(len(files_sample))
        for i in range(0,len(files_sample)):
            sp[i] = np.median(sample_image[256,256,i]) # This is for the first Bragg edge fitting
        #run once the fitting to check if everything works
        AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(myspectrum=sp, myrange=myrange, myTOF=mylambda_bin, est_pos=est_pos, est_sigma=est_sigma, est_alpha=est_alpha, bool_print=1, bool_average=0, bool_linear=0)

    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_heigth = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
#         print('processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                #print(i,j,' ciao')
                mysignal = np.zeros(len(files_sample))
                for ind in range(0,len(files_sample)):
                    mysignal[ind] = np.median(sample_image[i,j,ind])
                try:
                    edge_fit = AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(myspectrum=mysignal, myrange=myrange, myTOF=mylambda_bin, est_pos=est_pos, est_sigma=est_sigma, est_alpha=est_alpha, bool_print=0, bool_average=0, bool_linear=0)
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

    if(bool_print):
        plt.figure()
        plt.imshow(edge_position)
        plt.figure()
        plt.imshow(edge_width)
        plt.figure()
        plt.imshow(edge_heigth)        
    if(bool_save):
        np.save('edge_position.npy', edge_position)
        np.save('edge_height.npy', edge_heigth)
        np.save('edge_width.npy', edge_width)

    print("--- %s seconds ---" % (time.time() - start_time))
   
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width}
# dose is not finalized
# double check how the lambda range is picked
# 