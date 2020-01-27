
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import AdvancedBraggEdgeFitting_v2
# from PIL import Image

import tifffile
from tifffile import TiffFile

from skimage import io

import os, fnmatch
from os import listdir

from TOF_routines import tof2l
from TOF_routines import find_nearest

import time


#path to sample images
pathdata = '/media/carminati_c/Data2/LSP_Manuel/sample_binned_spotCleaned/'

#path to open beam images
pathob = '/media/carminati_c/Data2/LSP_Manuel/OB_binned_spotCleaned/'
files_sample = (sorted(fnmatch.filter(listdir(pathdata),'*.tif')))
files_ob = (sorted(fnmatch.filter(listdir(pathob),'*.tif')))

#load the spectrum
spectrum = np.loadtxt('/media/ws_niag/10_people/Morgano/RADEN_data_analysis/TEST6_000_Spectra.txt', usecols=0)


spectrum_binned = spectrum[0::18]
print(len(spectrum_binned))

#here che calibration parameters
t0=-0.0002618673892937752
L=18.961609065251505
mylambda_bin = tof2l(spectrum_binned,t0,L)

(pathob+files_ob[0])

ob = io.imread(pathob+files_ob[0])
np.shape(ob)

ob_image = np.zeros(np.shape(ob))
sample_image = np.zeros(np.shape(ob))


#read the images, in this example they are already selected in a neighborhood of the studied edge 
for i in range(0,155):
    ob_image[:,:,i] = io.imread(pathob+files_ob[i])
    sample_image[:,:,i] = io.imread(pathdata+files_sample[i])

#compute transmission image
trans_image = np.zeros([300,400,155])
for i in range(0,155):
    dose = np.median(ob_image[10:20,10:20,i])/np.median(sample_image[10:20,10:20,i]) # update the doseROI
    trans_image[:,:,i] = sample_image[:,:,i]/ob_image[:,:,i]*dose

sp=np.zeros(156)
for i in range(0,155):
    sp[i] = np.median(trans_image[44,204,i])


#here define the lambda range to be used for single edge fitting
myrange =[]
myrange.append(find_nearest(mylambda_bin, 3.7))
myrange.append(find_nearest(mylambda_bin, 4.4))
print(myrange)

# ------ Create the mask by thresholding ----
# mymask = trans_image[:,:,70]<0.5

# mymask[:40,:] = False
# mymask[:,380:] = False

# othermask = trans_image[:,:,70]<0.5
# othermask[:28,:] = False
# othermask[:,380:] = False

# little_mask = (othermask!=mymask)

# ---- Or load the mask -----
pathmask = '/pathtomask/'
filemask = 'mask.tif'
mymask = io.imread(pathmask+filemask)


est_sigma = 0.08
est_alpha = 0.08
est_pos = find_nearest(mylambda_bin[myrange[0]:myrange[1]], 4.05)
print(est_pos)
small_range = np.array([0, myrange[1]-myrange[0]-1])
small_lambda = mylambda_bin[myrange[0]:myrange[1]]
print(small_range)

#run once the fitting to check if everything works
AdvancedBraggEdgeFitting_v2.AdvancedBraggEdgeFitting(sp[myrange[0]:myrange[1]], small_range, small_lambda, est_pos, est_sigma, est_alpha, True, False, False, True)

edge_position = np.zeros(np.shape(mymask))
edge_width = np.zeros(np.shape(mymask))
print(np.shape(edge_position))

#loop for all pixel position, where the mask is equal to one
start_time = time.time()
for i in range(0, np.shape(mymask)[0]):

    print('processing row n. ', i, 'of', np.shape(mymask)[0])

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
                if (len(edge_fit['pos_extrema'])==2):
                    edge_width[i,j] = small_lambda[edge_fit['pos_extrema'][1]]-small_lambda[edge_fit['pos_extrema'][0]]
                else:
                    edge_width[i,j]=-2.0
            except:
                print("Unexpected error at :", i, j)
                edge_position[i,j]= -2.0
                edge_width[i,j]=-2.0


# plt.figure()
# plt.imshow(edge_position)
# save the results
np.save("edge_position.npy", edge_position)


print("--- %s seconds ---" % (time.time() - start_time))

