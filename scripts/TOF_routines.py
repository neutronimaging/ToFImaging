import numpy as np

h=6.62607004e-34 #Planck constant [m^2 kg / s]
m=1.674927471e-27 #Neutron mass [kg]

def tof2l(tof, lambda0, L):
    l=lambda0+h/m*(tof)/(L)/1e-10
    return l

def l2tof(l, t0, L):
    tof=t0+(l*1e-10)*(L)*m/h
    return tof

def binning (mysignal, newsize):
    binned_signal = np.zeros(newsize)
    bin_size = int(len(mysignal)/newsize)
    for i in range(0, newsize):
        bin_value = np.median(mysignal[i*bin_size:i*bin_size+bin_size])
        binned_signal[i]=bin_value
    return (binned_signal)
    
def binning_resolution (mysignal, spectrum, d_spectrum): #make it so it matches
    spectrum_range = spectrum[end]-spectrum[0]
    n_range = np.round(spectrum_range/d_spectrum)
    new_spectrum = np.arange(spectrum[0],spectrum[end]+spectrum_range/n_range,spectrum_range/n_range)
    if(len(np.shape(mysignal))==3):
        if(len(spectrum)!=np.shape(mysignal)[2]):
            print('Length of spectrum does not match signal size')
        binned_signal=np.zeros([np.shape(mysignal)[0], np.shape(mysignal)[1], len(new_spectrum)])
        for i in range(0,np.shape(mysignal)[0]):
            for j in range(0,np.shape(mysignal)[1]):
                binned_signal[i,j,:] = np.interp(new_spectrum,spectrum,mysignal[i,j,:])
                
    if(len(np.shape(mysignal))==2):
        if(len(spectrum)!=np.shape(mysignal)[1]):
            print('Length of spectrum does not match signal size')
        binned_signal=np.zeros([np.shape(mysignal)[0], len(new_spectrum)])
        for i in range(0,np.shape(mysignal)[0]):
                binned_signal[i,:] = np.interp(new_spectrum,spectrum,mysignal[i,:])
                
    if(len(np.shape(mysignal))==1):
        if(len(spectrum)!=np.shape(mysignal)[0]):
            print('Length of spectrum does not match signal size')
        binned_signal = np.interp(new_spectrum,spectrum,mysignal)
    
    return (binned_signal, new_spectrum)

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
    
    if(len(np.shape(mysignal))==3):
        outsignal = np.zeros((np.shape(mysignal)[0], np.shape(mysignal)[1], np.shape(mysignal)[2]))
        for i in range(0,np.shape(mysignal)[2]):
            outsignal[:,:,i] = scipy.signal.convolve2d(mysignal,K,'same')
    else:
        outsignal = scipy.signal.convolve2d(mysignal,K,'same')
    return outsignal   

def transmission_image (I,I0,dose_mask_path=0):
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
