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
import TOF_routines
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
    """ Run edge_fitting on provided 3D image, the third dimension is supposed to be lambda(A) or TOF 
    #INPUTS:
    #pathdata: path to sample transmission images
    #filemask: path to the .tif mask of the region to fit. if not specified runs over the whole area
    #pathspectrum: path to the .txt spectrum (in tof)
    #t0_cal: t0 for conversion from tof to lambda
    #L_cal: L for conversion from tof to lambda
    #lambda_range: lambda range for single edge fitting (in Angstrom): 2 values array
    #est_pos: expected bragg edge position (in A), if unset is found automatically
    #est_sigma: expected Gaussian broadening
    #est_alpha: expected decay constant (moderator property))
    #bool_save: flag for data saving
    #bool_print: flag for image printing
    
    ##OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    #'edge_position' : edge position
    #'edge_height': edge height
    #'edge_width': edge width
    """ 
    
    debug_flag = 1
    
    files_sample = (sorted(fnmatch.filter(listdir(pathdata),'*.tif')))
    files_ob = (sorted(fnmatch.filter(listdir(pathob),'*.tif')))
    #files_sample = (sorted(fnmatch.filter(listdir(pathdata),'*.fits')))
    #files_ob = (sorted(fnmatch.filter(listdir(pathob),'*.fits')))
    
    #load the spectrum
    spectrum = np.loadtxt(pathspectrum, usecols=0)
    spectrum_binned = spectrum[0::binning]
    
    mylambda_bin = tof2l(spectrum_binned,t0_cal,L_cal)

    ob = io.imread(pathob+'/'+files_ob[0])
    #ob = fits.getdata(pathob+'/'+files_ob[0])
    
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
    edge_height = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
        #print('processing row n. ', i, 'of', np.shape(mymask)[0])
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
                    edge_height[i,j] = edge_fit['height']
                    if (len(edge_fit['pos_extrema'])==2):
                        edge_width[i,j] = small_lambda[edge_fit['pos_extrema'][1]]-small_lambda[edge_fit['pos_extrema'][0]]
                    else:
                        edge_width[i,j] = -2.0
                except:
                    print("Unexpected error at :", i, j)
                    edge_position[i,j] = -2.0
                    edge_width[i,j] = -2.0
                    edge_height[i,j] = -2.0

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

    print("--- %s seconds ---" % (time.time() - start_time))
   
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width}

def image_edge_fitting_Ttof(pathdata, pathspectrum, lambda_range, filemask=0, t0_cal=-0.0002618673892937752, L_cal=18.961609065251505, est_pos=0, est_sigma=1, est_alpha=1, bool_save=0, bool_print=1):
    """ Run edge_fitting on provided 3D image, the third dimension is supposed to be lambda(A) or TOF 
    #INPUTS:
    #pathdata: path to sample transmission images
    #filemask: path to the .tif mask of the region to fit. if not specified runs over the whole area
    #pathspectrum: path to the .txt spectrum (in tof)
    #t0_cal: t0 for conversion from tof to lambda
    #L_cal: L for conversion from tof to lambda
    #lambda_range: lambda range for single edge fitting (in Angstrom): 2 values array
    #est_pos: expected bragg edge position (in A), if unset is found automatically
    #est_sigma: expected Gaussian broadening
    #est_alpha: expected decay constant (moderator property))
    #bool_save: flag for data saving
    #bool_print: flag for image printing
    
    ##OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    #'edge_position' : edge position
    #'edge_height': edge height
    #'edge_width': edge width
    """ 
   
    debug_flag = 0
    
    files_sample = (sorted(fnmatch.filter(listdir(pathdata),'*.tif')))
    
    #load the spectrum
    spectrum = np.loadtxt(pathspectrum, usecols=0)
    
    mylambda_bin = tof2l(spectrum,t0_cal,L_cal)

    ob = io.imread(pathdata+'/'+files_sample[0])
    
    sample_image = np.zeros([np.shape(ob)[0], np.shape(ob)[1], len(files_sample)])
    
    if(filemask):
        mymask = io.imread(filemask)
        if( [np.shape(sample_image)[0], np.shape(sample_image)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    else:
        mymask = np.ones([np.shape(sample_image)[0], np.shape(sample_image)[1]])

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
        AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(myspectrum=sp, myrange=myrange, myTOF=mylambda_bin, est_pos=est_pos, est_sigma=est_sigma, est_alpha=est_alpha, bool_print=debug_flag, bool_average=0, bool_linear=0)

    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
        #print('processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                print(i,j)
                mysignal = np.zeros(len(files_sample))
                for ind in range(0,len(files_sample)):
                    mysignal[ind] = np.median(sample_image[i,j,ind])
                try:
                    edge_fit = AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(myspectrum=mysignal, myrange=myrange, myTOF=mylambda_bin, est_pos=est_pos, est_sigma=est_sigma, est_alpha=est_alpha, bool_print=debug_flag, bool_average=0, bool_linear=0)
                    edge_position[i,j] = edge_fit['t0']
                    edge_height[i,j] = edge_fit['height']
                    if (len(edge_fit['pos_extrema'])==2):
                        edge_width[i,j] = mylambda_bin[edge_fit['pos_extrema'][1]]-mylambda_bin [edge_fit['pos_extrema'][0]]
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
   
    return {'edge_position' : edge_position, 'edge_height': edge_height, 'edge_width': edge_width}

def image_edge_fitting_Ttof_gauss(pathdata, pathspectrum, lambda_range, filemask=0, t0_cal=-0.0002618673892937752, L_cal=18.961609065251505, est_pos = 0, bool_save=0, bool_print=1):
    """ Run edge_fitting on provided 3D image, the third dimension is supposed to be lambda(A) or TOF     
    ##INPUTS:
    #pathdata: path to sample transmission images
    #filemask: path to the .tif mask of the region to fit. if not specified runs over the whole area
    #pathspectrum: path to the .txt spectrum (in tof)
    #t0_cal: t0 for conversion from tof to lambda
    #L_cal: L for conversion from tof to lambda
    #lambda_range: lambda range for single edge fitting (in Angstrom): 2 values array
    #est_pos: expected bragg edge position (in A), if unset is found automatically
    #bool_save: flag for data saving
    #bool_print: flag for image printing
    
    ##OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    #'edge_position' : edge position
    #'edge_height': edge height
    #'edge_width': edge width
    """ 
    
    debug_flag = 1
    
    files_sample = (sorted(fnmatch.filter(listdir(pathdata),'*.tif')))
    
    #load the spectrum
    spectrum = np.loadtxt(pathspectrum, usecols=0)
    
    mylambda_bin = tof2l(spectrum,t0_cal,L_cal)

    ob = io.imread(pathdata+'/'+files_sample[0])
    
    sample_image = np.zeros([np.shape(ob)[0], np.shape(ob)[1], len(files_sample)])
    
    if(filemask):
        mymask = io.imread(filemask)
        if( [np.shape(sample_image)[0], np.shape(sample_image)[1]] != [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    else:
        mymask = np.ones([np.shape(sample_image)[0], np.shape(sample_image)[1]])

    #read the images, in this example they are already selected in a neighborhood of the studied edge 
    for i in range(0,len(files_sample)):
        sample_image[:,:,i] = io.imread(pathdata+'/'+files_sample[i])/256

    #here define the lambda range to be used for single edge fitting
    myrange = []
    myrange.append(find_nearest(mylambda_bin, lambda_range[0])) # 3.7
    myrange.append(find_nearest(mylambda_bin, lambda_range[1])) # 4.2
    if(debug_flag): #testing on a single pixel    
        plt.imshow(sample_image.sum(axis=2))
        plt.show()
        plt.close()
        sp = np.zeros(len(files_sample))
        for i in range(0,len(files_sample)):
            sp[i] = np.median(sample_image[256,256,i]) # This is for the first Bragg edge fitting
        #run once the fitting to check if everything works
        AdvancedBraggEdgeFitting.GaussianBraggEdgeFitting(myspectrum=sp, myTOF=mylambda_bin, myrange=myrange, est_pos = est_pos, bool_smooth=1, smooth_w = 3, smooth_n = 0, bool_print=1)

    median_image = TOF_routines.medianimage(Ttof)
    edge_position = np.zeros(np.shape(mymask))
    edge_width = np.zeros(np.shape(mymask))
    edge_height = np.zeros(np.shape(mymask))
    #loop for all pixel position, where the mask is equal to one
    start_time = time.time()
    for i in range(0, np.shape(mymask)[0]):
        #print('processing row n. ', i, 'of', np.shape(mymask)[0])
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i,j]):
                print(i,j)
                mysignal = np.zeros(len(files_sample))
                for ind in range(0,len(files_sample)):
                    mysignal[ind] = np.median(sample_image[i,j,ind])
                try:
                    edge_fit = AdvancedBraggEdgeFitting.GaussianBraggEdgeFitting(myspectrum=mysignal, myTOF=mylambda_bin, myrange=myrange, est_pos = est_pos, bool_smooth=1, smooth_w = 3, smooth_n = 0, bool_print=debug_flag)
                    edge_position[i,j] = edge_fit['t0']
                    edge_height[i,j] = edge_fit['edge_height']
                    edge_width[i,j] = edge_fit['edge_width']
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
    
def image_edge_fitting_Tlambda(Ttof, spectrum_l, lambda_range, filemask=0, auto_mask = True, mask_thresh = [0.05, 0.95], est_pos=0, est_sigma=1, est_alpha=1, bool_save=True, bool_print=True, debug_flag=False, debug_idx=[160,200], bool_average=False, bool_linear=False):    
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
                    edge_fit = AdvancedBraggEdgeFitting.AdvancedBraggEdgeFitting(myspectrum=mysignal, myrange=myrange, myTOF=spectrum_l, est_pos=est_pos, est_sigma=est_sigma, est_alpha=est_alpha, bool_print=debug_flag, bool_average=bool_average, bool_linear=bool_linear)
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
    
def image_edge_fitting_Tlambda_gauss(Ttof, spectrum_l, lambda_range, filemask=0, auto_mask = True, mask_thresh = [0.05, 0.95], est_pos=0, est_sigma=1, est_alpha=1, bool_smooth=False, smooth_w = 3, smooth_n = 0, bool_save=True, bool_print=True, debug_flag=False, debug_idx=[160,200]):        
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
        sp = np.zeros(len(spectrum_l))
        for i in range(0,len(spectrum_l)):
            sp[i] = np.median(Ttof[debug_idx[0],debug_idx[1],i]) # This is for the first Bragg edge fitting
        #run once the fitting to check if everything works
        AdvancedBraggEdgeFitting.GaussianBraggEdgeFitting(myspectrum=sp, myTOF=spectrum_l, myrange=myrange, est_pos = est_pos, bool_smooth=bool_smooth, smooth_w = smooth_w, smooth_n = smooth_n, bool_print=debug_flag)

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
                mysignal = np.zeros(len(spectrum_l))
                for ind in range(0,len(spectrum_l)):
                    mysignal[ind] = np.median(Ttof[i,j,ind])
                try:
                    edge_fit = AdvancedBraggEdgeFitting.GaussianBraggEdgeFitting(myspectrum=mysignal, myTOF=spectrum_l, myrange=myrange, est_pos = est_pos, bool_smooth=bool_smooth, smooth_w = smooth_w, smooth_n = smooth_n, bool_print=debug_flag)
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