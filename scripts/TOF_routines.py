import numpy as np

h=6.62607004e-34 #Planck constant [m^2 kg / s]
m=1.674927471e-27 #Neutron mass [kg]

def tof2l(tof, L, lambda_0 = 0, tof_0 = 0):
    if(lambda_0):
        l=lambda_0+h/m*(tof)/(L)/1e-10
    if(tof_0):
        l=0.3956*(tof*1e6 +tof_0)/(L*100) #converts L to cm and tof0 must be in ns
    return l

def l2tof(l, t0, L):
    tof=t0+(l*1e-10)*(L)*m/h
    return tof
    
#these tools average stacks of image into a single one (axis =2)    
def averageimage(imgs):
    img=imgs.mean(axis=2)
    return img
    
def medianimage(imgs):
    img=np.median(imgs,axis=2)
    return img
    
def weightedaverageimage(imgs,size=5):
    import scipy.ndimage
    dims=imgs.shape
    w=np.zeros(imgs.shape)
    M=size**2
    for i in np.arange(dims[2]) :
        f=scipy.ndimage.uniform_filter(imgs[:,:,i], size=size)*M
        f2=scipy.ndimage.uniform_filter(imgs[:,:,i]**2, size=size)*M
        sigma=np.sqrt(1/(M-1)*(f2-(f**2)/M))
        
        w[:,:,i]=1.0/sigma
        
    wsum=w.sum(axis=2)
    for i in np.arange(dims[2]) :
        w[:,:,i]=w[:,:,i]/wsum
        
    imgs=w*imgs
    img=imgs.sum(axis=2)
    
    return img
    
def binning (mysignal, newsize):
    binned_signal = np.zeros(newsize)
    bin_size = int(len(mysignal)/newsize)
    for i in range(0, newsize):
        bin_value = np.median(mysignal[i*bin_size:i*bin_size+bin_size])
        binned_signal[i]=bin_value
    return (binned_signal)
    
def binning_resolution (mysignal, spectrum, d_spectrum):
    spectrum_range = spectrum[-1]-spectrum[0]
    n_range = np.round(spectrum_range/d_spectrum)
    new_spectrum = np.arange(spectrum[0],spectrum[-1]+spectrum_range/n_range,spectrum_range/n_range)
    if(len(np.shape(mysignal))==3):
        if(len(spectrum)!=np.shape(mysignal)[2]):
            print('WARNING: Length of spectrum does not match signal size')
        binned_signal=np.zeros([np.shape(mysignal)[0], np.shape(mysignal)[1], len(new_spectrum)])
        for i in range(0,np.shape(mysignal)[0]):
            for j in range(0,np.shape(mysignal)[1]):
                binned_signal[i,j,:] = np.interp(new_spectrum,spectrum,mysignal[i,j,:])
                
    if(len(np.shape(mysignal))==2):
        if(len(spectrum)!=np.shape(mysignal)[1]):
            print('WARNING: Length of spectrum does not match signal size')
        binned_signal=np.zeros([np.shape(mysignal)[0], len(new_spectrum)])
        for i in range(0,np.shape(mysignal)[0]):
                binned_signal[i,:] = np.interp(new_spectrum,spectrum,mysignal[i,:])
                
    if(len(np.shape(mysignal))==1):
        if(len(spectrum)!=np.shape(mysignal)[0]):
            print('WARNING: Length of spectrum does not match signal size')
        binned_signal = np.interp(new_spectrum,spectrum,mysignal)
    
    return (binned_signal, new_spectrum)
    
def load_fits (pathdata, cut_last=0, bool_wavg = False):
    from astropy.io import fits
    import os, fnmatch
    from os import listdir
    
    subfolders = sorted(listdir(pathdata))
    
    #load 1st subfolder and image to figure det shape and number of frames
    files = sorted(fnmatch.filter(listdir(pathdata+'\\'+subfolders[0]),'*.fits'))
    testim = fits.getdata(pathdata+'\\'+subfolders[0]+'\\'+files[0])
    det_shape = np.shape(testim)
    n_tof = len(files)
    data = np.zeros([det_shape[0], det_shape[1], n_tof-cut_last])
    
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
    data_t = np.zeros([det_shape[0], det_shape[1],len(subfolders)])
    for i in range(0,len(files)-cut_last):
        for j in range(0,len(subfolders)):
            files = sorted(fnmatch.filter(listdir(pathdata+'\\'+subfolders[j]),'*.fits'))
            data_t[:,:,j] = fits.getdata(pathdata+'\\'+subfolders[j]+'\\'+files[i])
        if (bool_wavg):
            data[:,:,i] = weightedaverageimage(data_t)    
        else:
            data[:,:,i] = medianimage(data_t)    
    return data

def transmission_normalization (I,I0,dose_mask_path=0):
    from skimage import io
    if(np.shape(I)!=np.shape(I0)):
        print('The size of I and I0 does not match. Check input data')
    if(len(np.shape(I))!=2):
        print('The input data is not an image')
    
    dose = 1
    if(dose_mask_path):
        dose_mask = io.imread(dose_mask_path)
        if(np.shape(I)!=np.shape(dose_mask)):
            print('The size of the dose mask does not match the signal. Check input data')
        dose_mask[dose_mask!=1]=np.nan
        dose = np.median(I0)/np.median(I)
    
    T = np.divide(I,I0)*dose
    T[T>1] = 1
    T[np.isinf(T)] = np.nan
    T[T<0] = np.nan
    return T

def interp_image_T (mydata):
    from scipy import interpolate
    mydata = np.double(mydata)
    mydata[mydata==0] = np.nan #can't be 0 transmission
    mydata[np.isinf(mydata)] = np.nan #caused by 0 counts in open beam
    
    x = np.arange(0,np.shape(mydata)[1])
    y = np.arange(0,np.shape(mydata)[0])
    
    mydata = np.ma.masked_invalid(mydata)
    xx, yy = np.meshgrid(x, y)
    x1 = xx[~mydata.mask]
    y1 = yy[~mydata.mask]
    newdata =  mydata[~mydata.mask]
    
    interp_data = interpolate.griddata((x1,y1), newdata.ravel(), (xx, yy), method='cubic')
    return interp_data
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return (idx)

def find_first(array, value):
    for idx in range(0, len(array)):
        if array[idx]>=value:
            break
    return idx

def find_last(array, value):
    for idx in range(0, len(array)):
        if idx>1:
            if ((array[idx]-array[idx-1])>=value+0.005):
                break
    return idx

def rotatedata(x,y,a):
    cosa=np.cos(a)
    sina=np.sin(a)
    x = x*cosa-y*sina
    y = x*sina+y*cosa
    return x,y

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
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
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def moving_average_1D (mysignal, kernel_size = 3, custom_kernel = 0):
    if(len(np.shape(mysignal))!=1):
        print('Data size is not 1D')
    if(custom_kernel):
        K = custom_kernel
    else:
        K = np.ones((kernel_size))
    K = K/np.sum(K)
    outsignal = np.convolve(mysignal,K,'same')
    return outsignal
    
def moving_average_2D (mysignal, kernel_size = 3, custom_kernel = 0):
    import scipy.signal
    if(len(np.shape(mysignal))!=3 | len(np.shape(mysignal))!=2):
        print('Data size is not either a 2D or ToF 2D')
    if(custom_kernel):
        K = custom_kernel
    else:
        K = np.ones((kernel_size,kernel_size))
    K = K/np.sum(K)
    
    if(len(np.shape(mysignal))==3): #if finds 3d matrix assume it's ToF data and apply to each tof frame
        outsignal = np.zeros((np.shape(mysignal)[0], np.shape(mysignal)[1], np.shape(mysignal)[2]))
        for i in range(0,np.shape(mysignal)[2]):
            outsignal[:,:,i] = scipy.signal.convolve2d(mysignal[:,:,i],K,'same')
    else:
        outsignal = scipy.signal.convolve2d(mysignal,K,'same')
    return outsignal   

def fullspectrum_T (path_sample, path_ob, cut_last=0):
    #load rawdata
    I = load_fits(path_sample,cut_last)
    I0 = load_fits(path_ob,cut_last)
    #normalize
    T = transmission_normalization(I.sum(axis=2),I0.sum(axis=2))
    #clean from nans/infs
    T = interp_image_T(T)   
    return(T)

def fullspectrum_im (path_data, cut_last=0):
    #load rawdata
    I = load_fits(path_data,cut_last)
    T = I.sum(axis=2)
    return(T)    
    
def load_routine (path_sample, path_ob, path_spectrum, cut_last=0, bin_size=0, d_spectrum = 0, dose_mask_path = 0, bool_lambda=False, L = 0, tof_0 = 0, lambda_0 = 0, bool_wavg = False, bool_interp = False):
    #load rawdata
    I = load_fits(path_sample,cut_last,bool_wavg)
    I0 = load_fits(path_ob,cut_last,bool_wavg)
    spectrum = np.loadtxt(path_spectrum, usecols=0)
    if(bool_lambda):
        spectrum = tof2l(spectrum, L, lambda_0, tof_0)
    #rebinning
    if(d_spectrum):
        I = binning_resolution(I,spectrum,d_spectrum)
        I0 = binning_resolution(I0,spectrum,d_spectrum)
    if(bin_size):
        I = binning(I,bin_size)
        I0 = binning(I0,bin_size)
    #normalize
    T = transmission_normalization(I,I0,dose_mask_path)
    #clean from nans/infs
    if(bool_interp):
        for i in range(0,np.shape(T)[2]):
            T[:,:,i] = interp_image_T(T[:,:,i])
        
    return{'T':T, 'spectrum':spectrum}
