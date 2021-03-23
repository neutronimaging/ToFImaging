import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

h = 6.62607004e-34  #Planck constant [m^2 kg / s]
m = 1.674927471e-27  #Neutron mass [kg]


#Unit Conversions
def Ang2MeV(Angstrom):
    return 81.82 / (Angstrom**2)


def MeV2Ang(MeV):
    return np.sqrt(81.82 / MeV)


#Calibration functions
def tof2l(tof, L, lambda_0=0, tof_0=0):
    if (lambda_0):
        l = lambda_0 + h / m * (tof) / (L) / 1e-10
    if (tof_0):
        l = 0.3956 * (tof * 1e6 + tof_0) / (
            L * 100)  #converts L to cm and tof0 must be in ns
    return l


def l2tof(l, t0, L):
    tof = t0 + (l * 1e-10) * (L) * m / h
    return tof


def tof2l_t0k(t, t0, k):
    return t0 + k * t


#Multiple frame merging tools, useful when an acquisition is split into sub-acquisitions (axis=2)
def averageimage(imgs):
    img = imgs.mean(axis=2)
    return img


def medianimage(imgs):
    img = np.median(imgs, axis=2)
    return img


def weightedaverageimage(imgs, size=5):
    import scipy.ndimage
    dims = imgs.shape
    w = np.zeros(imgs.shape)
    M = size**2
    for i in np.arange(dims[2]):
        f = scipy.ndimage.uniform_filter(imgs[:, :, i], size=size) * M
        f2 = scipy.ndimage.uniform_filter(imgs[:, :, i]**2, size=size) * M
        sigma = np.sqrt(1 / (M - 1) * (f2 - (f**2) / M))

        w[:, :, i] = 1.0 / sigma

    wsum = w.sum(axis=2)
    for i in np.arange(dims[2]):
        w[:, :, i] = w[:, :, i] / wsum

    imgs = w * imgs
    img = imgs.sum(axis=2)

    return img


#Rebinning/averaging tools
def binning_ndarray(mysignal, newsize):
    #Remnant from Chiara's code: Rebins an ndarray into the newsize.
    binned_signal = np.zeros(newsize)
    bin_size = int(len(mysignal) / newsize)
    for i in range(0, newsize):
        bin_value = np.median(mysignal[i * bin_size:i * bin_size + bin_size])
        binned_signal[i] = bin_value
    return (binned_signal)


def spectral_binning_resolution(mysignal, spectrum, d_spectrum):
    # Rebins an (spectrum) or (x,spectrum) or (x,y,spectrum) matrix to match a new bin_width (e.g. instrumental resolution)
    spectrum_range = spectrum[-1] - spectrum[0]
    n_range = np.round(spectrum_range / d_spectrum)
    new_spectrum = np.arange(spectrum[0],
                             spectrum[-1] + spectrum_range / n_range,
                             spectrum_range / n_range)
    if (len(np.shape(mysignal)) == 3):
        if (len(spectrum) != np.shape(mysignal)[2]):
            print('WARNING: Length of spectrum does not match signal size')
        binned_signal = np.zeros(
            [np.shape(mysignal)[0],
             np.shape(mysignal)[1],
             len(new_spectrum)])
        for i in range(0, np.shape(mysignal)[0]):
            for j in range(0, np.shape(mysignal)[1]):
                binned_signal[i, j, :] = np.interp(new_spectrum, spectrum,
                                                   mysignal[i, j, :])

    if (len(np.shape(mysignal)) == 2):
        if (len(spectrum) != np.shape(mysignal)[1]):
            print('WARNING: Length of spectrum does not match signal size')
        binned_signal = np.zeros([np.shape(mysignal)[0], len(new_spectrum)])
        for i in range(0, np.shape(mysignal)[0]):
            binned_signal[i, :] = np.interp(new_spectrum, spectrum,
                                            mysignal[i, :])

    if (len(np.shape(mysignal)) == 1):
        if (len(spectrum) != np.shape(mysignal)[0]):
            print('WARNING: Length of spectrum does not match signal size')
        binned_signal = np.interp(new_spectrum, spectrum, mysignal)

    return (binned_signal, new_spectrum)


def moving_average_1D(mysignal, kernel_size=3, custom_kernel=np.ndarray([0])):
    # Moving average by kernel convolution to a ndarray
    if (len(np.shape(mysignal)) != 1):
        print('Data size is not 1D')
    if (custom_kernel.any()):
        K = custom_kernel
    else:
        K = np.ones((kernel_size))
    K = K / np.sum(K)
    outsignal = np.convolve(mysignal, K, 'same')
    return outsignal


def moving_average_2D(mysignal, box_kernel=[], custom_kernel=np.ndarray([0])):
    # Moving average by kernel convolution to an image (2D)
    # !! If it finds 3d matrix assume it's ToF data and apply to each tof frame !!
    import scipy.signal

    if (len(np.shape(mysignal)) != 3 | len(np.shape(mysignal)) != 2):
        print('Data size is not either a 2D or ToF 2D')
    if (custom_kernel.any()):
        K = custom_kernel
    elif (any(box_kernel)):
        M = np.max(box_kernel)
        m = np.min(box_kernel)
        d = np.argmin(box_kernel)
        K = np.zeros((M, M))
        if (M % 2 == 0):
            if (m % 2 == 0):
                K[np.int(M / 2 - m / 2):np.int(M / 2 + m / 2), :] = 1
            else:
                K[np.int(M / 2 -
                         np.floor(m / 2)):np.int(M / 2 +
                                                 np.floor(m / 2)), :] = 1
                K[np.int(M / 2 - np.floor(m / 2)) - 1, :] = 0.5
                K[np.int(M / 2 + np.floor(m / 2)), :] = 0.5
        else:
            if (m % 2 == 0):
                K[np.int(M / 2 - (m - 1) / 2):np.int(M / 2 +
                                                     (m - 1) / 2), :] = 1
                K[np.int(M / 2 - (m - 1) / 2) - 1, :] = 0.5
                K[np.int(M / 2 + (m - 1) / 2), :] = 0.5
            else:
                K[np.int(M / 2 - m / 2):np.int(M / 2 + m / 2), :] = 1
        if (d == 1):
            K = np.transpose(K)
    K = K / np.sum(K)

    if (
            len(np.shape(mysignal)) == 3
    ):  #if finds 3d matrix assume it's ToF data and apply to each tof frame
        outsignal = np.zeros((np.shape(mysignal)[0], np.shape(mysignal)[1],
                              np.shape(mysignal)[2]))
        for i in tqdm(range(0, np.shape(mysignal)[2])):
            outsignal[:, :, i] = scipy.signal.convolve2d(
                mysignal[:, :, i], K, 'same')
    else:
        outsignal = scipy.signal.convolve2d(mysignal, K, 'same')
    return outsignal


def DataFiltering(mysignal, BoxKernel=[], GaussianKernel=[], bool_print=False):
    import scipy.ndimage
    # ERRORS in dataset
    if (len(np.shape(mysignal)) == 1):
        print('Data must be 2D image or 3D volume, or TOF stack of 2D images.')
        return
    if ((len(BoxKernel) == 1 or len(GaussianKernel) == 1)):
        print('Kernels must be at least of size 2.')
        return
    if (len(np.shape(mysignal)) == 2
            and (len(BoxKernel) == 3 or len(GaussianKernel) == 3)):
        print('Data is 2d but filtering kernel is 3D.')
        return

    if (bool_print):
        plt.figure()
        plt.subplot(1, 2, 1),
        if (len(np.shape(mysignal)) == 3):
            plt.imshow(np.nanmean(
                mysignal, axis=2)), plt.title('Input image'), plt.colorbar()
        else:
            plt.imshow(mysignal), plt.title('Input image'), plt.colorbar()

    # TOF data (3D), 2D kernel
    if (len(np.shape(mysignal)) == 3
            and (len(BoxKernel) == 2 or len(GaussianKernel) == 2)):
        print(
            'Data is 3D but filtering kernel is 2D. Applying filter to each slice of the data (third dimension).'
        )
        outsignal = np.zeros((np.shape(mysignal)[0], np.shape(mysignal)[1],
                              np.shape(mysignal)[2]))
        if (any(BoxKernel)):
            kernel = np.ones((BoxKernel[0], BoxKernel[1]))
            kernel = kernel / np.sum(np.ravel(kernel))
            for i in tqdm(range(0, np.shape(mysignal)[2])):
                outsignal[:, :, i] = scipy.ndimage.convolve(
                    mysignal[:, :, i], kernel)

            if (bool_print):
                plt.subplot(1, 2, 2),
                plt.imshow(np.nanmean(outsignal,
                                      axis=2)), plt.title('Output image')
                plt.colorbar(), plt.tight_layout(), plt.show(), plt.close()
            return outsignal
        if (any(GaussianKernel)):
            for i in tqdm(range(0, np.shape(mysignal)[2])):
                outsignal[:, :, i] = scipy.ndimage.gaussian_filter(
                    mysignal[:, :, i], GaussianKernel)

            if (bool_print):
                plt.subplot(1, 2, 2),
                plt.imshow(np.nanmean(outsignal,
                                      axis=2)), plt.title('Output image')
                plt.colorbar(), plt.tight_layout(), plt.show(), plt.close()
            return outsignal

    # TOF data (3D), 3D kernel
    if (len(np.shape(mysignal)) == 3
            and (len(BoxKernel) == 3 or len(GaussianKernel) == 3)):
        print(
            'Data and filtering kernel are 3D. Applying 3D filter convolution.'
        )
        if (any(BoxKernel)):
            kernel = np.ones((BoxKernel[0], BoxKernel[1], BoxKernel[2]))
            kernel = kernel / np.sum(np.ravel(kernel))
            outsignal = scipy.ndimage.convolve(mysignal, kernel)
            if (bool_print):
                plt.subplot(1, 2, 2),
                plt.imshow(np.nanmean(outsignal,
                                      axis=2)), plt.title('Output image')
                plt.colorbar(), plt.tight_layout(), plt.show(), plt.close()
            return outsignal
        if (any(GaussianKernel)):
            outsignal = scipy.ndimage.gaussian_filter(mysignal, GaussianKernel)
            if (bool_print):
                plt.subplot(1, 2, 2),
                plt.imshow(np.nanmean(outsignal,
                                      axis=2)), plt.title('Output image')
                plt.colorbar(), plt.tight_layout(), plt.show(), plt.close()
            return outsignal

    # image data (2D), 2D kernel
    if (len(np.shape(mysignal)) == 2
            and (len(BoxKernel) == 2 or len(GaussianKernel) == 2)):
        print(
            'Data and filtering kernel are 2D. Applying 2D filter convolution.'
        )
        if (any(BoxKernel)):
            kernel = np.ones((BoxKernel[0], BoxKernel[1]))
            kernel = kernel / np.sum(np.ravel(kernel))
            outsignal = scipy.ndimage.convolve(mysignal, kernel)
            if (bool_print):
                plt.subplot(1, 2, 2),
                plt.imshow(outsignal), plt.title('Output image')
                plt.tight_layout(), plt.show(), plt.close()
            return outsignal
        if (any(GaussianKernel)):
            outsignal = scipy.ndimage.gaussian_filter(mysignal, GaussianKernel)
            if (bool_print):
                plt.subplot(1, 2, 2),
                plt.imshow(outsignal), plt.title('Output image')
                plt.tight_layout(), plt.show(), plt.close()
            return outsignal

    return


def spatial_image_rebinning(image, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and 
        new axes must divide old ones.

    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)

    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]

    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if image.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(
            image.shape, new_shape))
    compression_pairs = [(d, c // d) for d, c in zip(new_shape, image.shape)]
    flattened = [l for p in compression_pairs for l in p]
    image = image.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(image, operation)
        image = op(-1 * (i + 1))
    return image


def spatial_discrete_rebinning(image, rebinning_order=2, operation='sum'):
    def rebin(a, new_shape, operation='mean'):
        M, N = a.shape
        m = new_shape[0]
        n = new_shape[1]
        if (operation == 'sum'):
            a = a.reshape((m, np.int(M / m), n, np.int(N / n))).sum(3).sum(1)
        if (operation == 'mean'):
            a = a.reshape((m, np.int(M / m), n, np.int(N / n))).mean(3).mean(1)
        return a

    f_dim = (np.int(np.shape(image)[0] / rebinning_order),
             np.int(np.shape(image)[0] / rebinning_order))
    if ((np.shape(image)[0] % rebinning_order) != 0):
        print(
            'WARNING: the rebinning order does not divide the image dimension, please use a different function.'
        )

    if (
            len(np.shape(image)) == 3
    ):  #if finds 3d matrix assume it's ToF data and apply to each tof frame
        outsignal = np.zeros((f_dim[0], f_dim[1], np.int(np.shape(image)[2])))
        for i in tqdm(range(0, np.shape(image)[2])):
            outsignal[:, :, i] = rebin(image[:, :, i], f_dim, operation)
    else:
        outsignal = rebin(image, f_dim, operation)
    return outsignal


def tof_image_rebinning(image, spectrum, rebinning_order, operation='mean'):
    tof_n = np.shape(image)[2]
    tof_n_new = np.int(np.round(tof_n / rebinning_order))
    image_out = np.zeros((np.shape(image)[0], np.shape(image)[1], tof_n_new))
    spectrum_out = np.zeros((tof_n_new))
    for i in tqdm(range(0, tof_n_new)):
        if (operation == 'sum'):
            image_out[:, :, i] = np.nansum(
                image[:, :, rebinning_order * i:rebinning_order * i +
                      rebinning_order],
                axis=2)
            spectrum_out[i] = np.nansum(
                spectrum[rebinning_order * i:rebinning_order * i +
                         rebinning_order])
        elif (operation == 'mean'):
            image_out[:, :, i] = np.nanmean(
                image[:, :, rebinning_order * i:rebinning_order * i +
                      rebinning_order],
                axis=2)
            spectrum_out[i] = np.nanmean(
                spectrum[rebinning_order * i:rebinning_order * i +
                         rebinning_order])
    return (image_out, spectrum_out)


def savitzky_golay(y, window_size=5, order=1, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range]
                for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


#Normalization tools
def transmission_normalization(I, I0, mask_dose=np.ndarray([0])):
    #Normalization to transmission, with dose correction if mask_dose is given. negative transmission and negative attenuation (T>1) are replaced with NaNs.
    from skimage import io
    if (np.shape(I) != np.shape(I0)):
        print('The size of I and I0 does not match. Check input data')
    if (len(np.shape(I)) != 2):
        print(
            'The input data is not an image. Assuming is a TOF stack of images.'
        )

    dose = 1
    if (mask_dose.any()):
        dose = np.median(np.multiply(I0, mask_dose), axis=(1, 0)) / np.median(
            np.multiply(I, mask_dose), axis=(1, 0))

    T = np.divide(I * dose, I0)
    T[T > 1] = 1
    T[np.isinf(T)] = np.nan
    T[T < 0] = np.nan
    return T


#Index seeking functions
def find_nearest(array, value):
    #finds the nearest index to the value in the array
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return (idx)


def find_first(array, value):
    #finds the first index where the array is higher than the value from the start
    for idx in range(0, len(array)):
        if array[idx] >= value:
            break
    return idx


def find_last(array, value):
    #finds the last index where the array is higher than the value from the start
    for idx in range(0, len(array)):
        if idx > 1:
            if ((array[idx] - array[idx - 1]) >= value + 0.005):
                break
    return idx


def rotatedata(x, y, a):
    #rotates x,y by a(rad)
    cosa = np.cos(a)
    sina = np.sin(a)
    x = x * cosa - y * sina
    y = x * sina + y * cosa
    return x, y


#Quick image plotting functions
def fullspectrum_T(path_sample, path_ob, cut_last=0):
    #load rawdata
    I = load_fits(path_sample, cut_last)
    I0 = load_fits(path_ob, cut_last)
    #normalize
    T = transmission_normalization(I.sum(axis=2), I0.sum(axis=2))
    return (T)


def fullspectrum_im(path_data, cut_last=0):
    #load rawdata
    I = load_fits(path_data, cut_last)
    T = I.sum(axis=2)
    return (T)


#Segmentation tools
def SpectralSegmentation(T_tof,
                         clusters,
                         spectrum=[],
                         spectrum_range=[],
                         bool_print=1):
    from sklearn.cluster import KMeans

    if (spectrum_range):
        idx_low = find_nearest(spectrum, spectrum_range[0])
        idx_high = find_nearest(spectrum, spectrum_range[1])
        T_tof = T_tof[:, :, idx_low:idx_high]
        spectrum = spectrum[idx_low:idx_high]

    Tarray = np.reshape(
        T_tof, [np.shape(T_tof)[0] * np.shape(T_tof)[1],
                np.shape(T_tof)[2]])
    Tarray[np.isnan(Tarray)] = 0
    Tarray[np.isinf(Tarray)] = 0
    kmeans = KMeans(n_clusters=clusters, random_state=0).fit(Tarray)

    T_segmented = np.reshape(
        kmeans.labels_,
        [np.shape(T_tof)[0], np.shape(T_tof)[1]])
    spectra = kmeans.cluster_centers_

    if (bool_print):
        plt.imshow(T_segmented)
        plt.title('Segmented image')
        plt.colorbar()
        plt.show()
        plt.close()
        if (len(spectrum) == 0):
            plt.plot(np.transpose(spectra))
        else:
            plt.plot(spectrum, np.transpose(spectra))
        plt.title('Segmented spectra')
        plt.show()
        plt.close()

    return {'T_segmented': T_segmented, 'spectra': spectra}


#Loading routines
def load_fits(pathdata, cut_last=0, bool_wavg=False):
    # Load stacks of TOF into 3D matrix (x,y,TOF). Requires subfolders in the pathdata with a subfolder for each repetition that is merged.
    # It is often the case that a long acquisition is split into multiple shorter acquisition to be merged (e.g. 4 hours = 4x1 hour)
    from astropy.io import fits
    import os, fnmatch
    from os import listdir

    subfolders = sorted(listdir(pathdata))

    #load 1st subfolder and image to figure det shape and number of frames
    files = sorted(
        fnmatch.filter(listdir(pathdata + '\\' + subfolders[0]), '*.fits'))
    testim = fits.getdata(pathdata + '\\' + subfolders[0] + '\\' + files[0])
    det_shape = np.shape(testim)
    n_tof = len(files)
    data = np.zeros([det_shape[0], det_shape[1], n_tof - cut_last])

    # load all
    # for j in range(0,len(subfolders)):
    # files = sorted(fnmatch.filter(listdir(pathdata+'\\'+subfolders[j]),'*.fits'))
    # if(len(files)!=n_tof):
    # print('WARNING: Data has non-equal tof frames!')
    # for i in range(0,len(files)-cut_last):
    # data_t = fits.getdata(pathdata+'\\'+subfolders[j]+'\\'+files[i])
    # if(np.shape(data_t)!=det_shape):
    # print('WARNING: Data has different Detector shape!')
    # data[:,:,i] += data_t

    # would maybe nice to add check before this stage that the signals are similar before merging
    # data /= len(subfolders)
    data_t = np.zeros([det_shape[0], det_shape[1], len(subfolders)])
    print('Loading ' + pathdata)
    for i in tqdm(range(0, len(files) - cut_last)):
        for j in range(0, len(subfolders)):
            files = sorted(
                fnmatch.filter(listdir(pathdata + '\\' + subfolders[j]),
                               '*.fits'))
            data_t[:, :, j] = fits.getdata(pathdata + '\\' + subfolders[j] +
                                           '\\' + files[i])
        if (bool_wavg):
            data[:, :, i] = weightedaverageimage(data_t)
        else:
            data[:, :, i] = medianimage(data_t)
    return data


def load_routine(path_sample,
                 path_ob,
                 path_spectrum,
                 cut_last=0,
                 dose_mask=np.ndarray([0]),
                 bool_lambda=False,
                 L=0,
                 tof_0=0,
                 lambda_0=0,
                 bool_wavg=False,
                 tofrebin_order=0):
    # Full loading routine. Load sample and open beam and normalize to TOF transmission T=(x,y,TOF). The tof spectrum is loaded as well and converted to lambda, when asked.

    #load rawdata
    I = load_fits(path_sample, cut_last, bool_wavg)
    I0 = load_fits(path_ob, cut_last, bool_wavg)
    spectrum = np.loadtxt(path_spectrum, usecols=0)
    if (bool_lambda):
        spectrum = tof2l(spectrum, L, lambda_0, tof_0)
    #rebinning
    if (tofrebin_order):
        sp = spectrum
        I = tof_image_rebinning(I, sp, tofrebin_order)
        I0 = tof_image_rebinning(I0, sp, tofrebin_order)
        spectrum = np.zeros((np.shape(I)[2]))
        for i in range(0, np.shape(I)[2]):
            spectrum[i] = sp[i * tofrebin_order]

    #normalize
    T = transmission_normalization(I, I0, dose_mask)

    return {'T': T, 'spectrum': spectrum}