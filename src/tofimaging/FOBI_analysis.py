import numpy as np
import tofimaging.ReductionTools as rt
from tqdm import tqdm

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

def chopper_time_delays_generator(time,
                                  nslits=8,
                                  nrep=2,
                                  mode='pseudorandom',
                                  rng=0.25):
    """
    Generates time delays for a multi-slits chopper:
    
    INPUTS:
    time = time-of-flight bins 
    nslits = number of slits per pattern repetitions
    nrep = number of slit pattern repetitions
    mode = how slits are generated:
            - even = evenly distributed
            - random = randomly distributed
            - pseudorandom = for each pattern the slits are evenly distributed and then offset randomly scaled by the "rng" parameter
    
    OUTPUT:
    D =  array with discretized dirac deltas at the chopper time delays.
    """
    Nt = int(np.shape(time)[0])
    angles = np.zeros(nslits)

    if mode == 'even':
        angles = np.linspace(0, 360, nslits + 1) / nrep
        angles = angles[1:-2]
        ang_t = angles
    if mode == 'random':
        angles = 360 * np.random.normal(size=nslits) / nrep
        ang_t = angles
    if mode == 'pseudorandom':
        for i in range(0, nslits):
            angles[i] = i * 360 / (nslits * nrep) + np.random.normal(
            ) * 360 / (nslits * nrep) * rng
            ang_t = angles

    if nrep > 1:
        for rep in range(0, nrep - 1):
            angles = np.append(angles, ang_t + (rep + 1) * 360 / nrep)

    shifts = Nt * angles / 360
    D = np.zeros((Nt, 1))
    for i in range(0, int(np.shape(shifts)[0])):
        sfloor = np.floor(shifts[i])
        rest = shifts[i] - sfloor
        sfloor = int(sfloor + 1)
        D[sfloor] = 1 - rest
        D[sfloor + 1] = rest
    return D


def time_delays_poldi(time):
    """
    Generates time delays for the POLDI chopper:
    
    INPUTS:
    time = time-of-flight bins
    
    OUTPUT:
    D =  array with discretized dirac deltas at the chopper time delays.
    """
    Nt = int(np.shape(time)[0])
    angles = np.array(
        [0, 9.363, 21.475, 37.039, 50.417, 56.664, 67.422, 75.406])
    angles = 90 - angles
    angles = np.flipud(angles)
    angles = angles - angles[0]
    # angles = a

    angles = np.concatenate((angles, angles + 1.0 * 90.0, angles + 2.0 * 90.0,
                             angles + 3.0 * 90.0)) / 360.0
    shifts = Nt * angles
    D = np.zeros((Nt, 1))
    for i in range(0, np.shape(shifts)[0]):
        sfloor = np.floor(shifts[i])
        rest = shifts[i] - sfloor
        sfloor = np.int(sfloor + 1)
        D[sfloor] = 1 - rest
        D[sfloor + 1] = rest

    D = np.squeeze(D)
    return D


def time_delays_4x10(time):
    """
    Generates time delays for the 4x10 disk chopper:
    
    INPUTS:
    time = time-of-flight bins
    
    OUTPUT:
    D =  array with discretized dirac deltas at the chopper time delays.
    """
    Nt = int(np.shape(time)[0])
    angles = np.array([
        7.579, 13.316, 18.1825, 30.831, 36.201, 53.363, 61.1375, 67.332,
        76.674, 87.573
    ])
    angles = angles - angles[0]

    angles = np.concatenate((angles, angles + 1.0 * 90.0, angles + 2.0 * 90.0,
                             angles + 3.0 * 90.0)) / 360.0
    shifts = Nt * angles
    D = np.zeros((Nt, 1))
    for i in range(0, np.shape(shifts)[0]):
        sfloor = np.floor(shifts[i])
        rest = shifts[i] - sfloor
        sfloor = np.int(sfloor + 1)
        D[sfloor] = 1 - rest
        D[sfloor + 1] = rest

    D = np.squeeze(D)
    return D


def time_delays_5x8(time):
    """
    Generates time delays for the 4x10 disk chopper:
    
    INPUTS:
    time = time-of-flight bins
    
    OUTPUT:
    D =  array with discretized dirac deltas at the chopper time delays.
    """
    Nt = int(np.shape(time)[0])
    angles = np.array([4.81, 17.09, 24.47, 31.35, 36.48, 48.21, 59.64, 66.64])
    angles = angles - angles[0]

    angles = np.concatenate((angles, angles + 1.0 * 72.0, angles + 2.0 * 72.0,
                             angles + 3.0 * 72.0, angles + 4.0 * 72.0)) / 360.0
    shifts = Nt * angles
    D = np.zeros((Nt, 1))
    for i in range(0, np.shape(shifts)[0]):
        sfloor = np.floor(shifts[i])
        rest = shifts[i] - sfloor
        sfloor = np.int(sfloor + 1)
        D[sfloor] = 1 - rest
        D[sfloor + 1] = rest

    D = np.squeeze(D)
    return D


def time_delays_3x14(time):
    """
    Generates time delays for the 4x10 disk chopper:
    
    INPUTS:
    time = time-of-flight bins
    
    OUTPUT:
    D =  array with discretized dirac deltas at the chopper time delays.
    """
    Nt = int(np.shape(time)[0])
    angles = np.array([0.8958 10.4775 21.1331 32.2086 37.0404 43.7625 59.2190 65.4966 75.5918 85.4641 91.0146 98.9699 108.2814 113.3630])
    angles = angles - angles[0]
    angles = angles/120
    shifts = Nt * angles
    D = np.zeros((Nt, 1))
    for i in range(0, np.shape(shifts)[0]):
        sfloor = np.floor(shifts[i])
        rest = shifts[i] - sfloor
        sfloor = np.int(sfloor + 1)
        D[sfloor] = 1 - rest
        D[sfloor + 1] = rest

    D = np.squeeze(D)
    return D


def wiener_decorrelation(f, g, c=1e-1):
    """
    Perform the decorrelation of FOBI overlap patterns
    
    INPUTS:
    f =  measured signal
    g = time delays
    c = Wiener coefficient
    
    OUTPUT:
    H = reconstructed signal
    """
    F = np.fft.fft(f)
    G = np.fft.fft(g)
    arg = np.divide(F * G, (np.power(np.abs(G), 2) + c))
    H = np.fft.ifft(arg)
    return H


def wiener_deconvolution(f, g, c=1e-1):
    """
    Perform the decorrelation of FOBI overlap patterns
    
    INPUTS:
    f =  measured signal
    g = time delays
    c = Wiener coefficient
    
    OUTPUT:
    H = reconstructed signal
    """
    F = np.fft.fft(f)
    G = np.fft.fft(g)
    arg = np.divide(F * np.conj(G), (np.power(np.abs(G), 2) + c))
    H = np.fft.ifft(arg)
    return H


def interpolate_noreadoutgaps(y, t, tmax, nrep, plot_flag):
    # Reformat arrays
    y = np.squeeze(y)
    if y.shape[1] > 1:
        y = y.T
    t = np.squeeze(t)
    if t.shape[1] > 1:
        t = t.T

    # Get tof bin width, figure length, and merge overlaps
    replen = int(np.ceil(len(y) / nrep))
    t_tot = np.linspace(t[0], tmax, nrep * replen)
    y_int = np.interp(t_tot, t, y)
    y_overlap = np.zeros((replen, nrep))

    for i in range(nrep):
        y_overlap[:, i] = y_int[replen * i:replen * (i + 1)]
    
    y_merged = np.nanmean(y_overlap, axis=1)
    t_merged = t_tot[:replen]

    if plot_flag:
        plt.figure()
        for i in range(nrep):
            plt.plot(t[:replen], y_overlap[:, i])
        plt.plot(t[:replen], y_merged, 'k')
        plt.show()
    
    return y_merged, t_merged


def full_fobi_reduction(y,
                        y0,
                        t,
                        tmax,
                        nrep,
                        chopper_id,
                        c=1e-1,
                        roll_value=201,
                        bool_smooth=True,
                        SG_w=5,
                        SG_o=1):
    """ Performs FOBI reduction from open beam and sample spectrum    
    
    INPUTS:
    y = sample spectrum (I)
    y0 = open beam spectrum (I0)
    t = time-of-flight bins
    tmax = maximum time of flight (this parameter is dependent on the chopper frequency: tmax = 1/f  with f = chopper frequency)
    nrep = number of times the pattern is repeated (chopper)
    c = Wiener constant
    bool_roll = if this is activated the spectra are offset to have the minimum of the I0 spectrum as first time bin
    use for single spectra, but not for 2D spectra as it my yield an offset in the results
    SG_w = Savitzky Golay filter window
    SG_o = Savitzky Golay filter order

    OUTPUTS:
    y_fobi: FOBI reduced I
    y0_fobi: FOBI reduced I0
    T_fobi: FOBI reduced I/I0
    t_fobi: FOBI reduced TOF bins
    """
    if (bool_smooth):
        y0 = rt.savitzky_golay(y0, SG_w, SG_o)
        y = rt.savitzky_golay(y, SG_w, SG_o)

    [y, tn] = interp_noreadoutgaps(y, t, tmax, nrep)
    [y0, tn] = interp_noreadoutgaps(y0, t, tmax, nrep)
    if (chopper_id == 'poldi'):
        D = time_delays_poldi(tn)
    if (chopper_id == '4x10'):
        D = time_delays_4x10(tn)
    if (chopper_id == '5x8'):
        D = time_delays_5x8(tn)
    if (chopper_id == '3x14'):
        D = time_delays_3x14(tn)

    y0rec = wiener_deconvolution(y0, D, c)
    yrec = wiener_deconvolution(y, D, c)
    Trec = wiener_deconvolution(y/y0, D, c)

    if (roll_value == True):
        y0rec = np.roll(y0rec, -roll_value)
        yrec = np.roll(yrec, -roll_value)
        Trec = np.roll(Trec, -roll_value)

    t_fobi = tn

    return y0rec, yrec, Trec, t_fobi


def fobi_2d(I,
            I0,
            t,
            tmax,
            nrep,
            chopper_id,
            c=1e-1,
            roll_value=201,
            bool_smooth=True,
            SG_w=5,
            SG_o=1):
    """ Performs FOBI reduction from open beam and sample spectrum to a stack of images  
    
    INPUTS:
    y = sample spectrum (I)
    y0 = open beam spectrum (I0)
    t = time-of-flight bins
    tmax = maximum time of flight (this parameter is dependent on the chopper frequency: tmax = 1/f  with f = chopper frequency)
    nrep = number of times the pattern is repeated (chopper)
    c = Wiener constant
    roll_value = shift the retrieved fobi spectra by this value 
    (this is necessary because the output from Wiener decorrelation places the peak of 
    the spectrum approximately in the middle of the array) This depends on time of flight bins and type of choppper
    bool_smooth = Set to true to apply S-G filter
    SG_w = Savitzky Golay filter window
    SG_o = Savitzky Golay filter order

    OUTPUTS:
    y_fobi: FOBI reduced I
    y0_fobi: FOBI reduced I0
    T_fobi: FOBI reduced I/I0
    t_fobi: FOBI reduced TOF bins
    """
    #get reference
    yref = np.nanmean(np.nanmean(I, axis=0), axis=0)
    y0ref = np.nanmean(np.nanmean(I0, axis=0), axis=0)
    [y0ref, yref, Tref, tref] = full_fobi_reduction(yref,
                                                    y0ref,
                                                    t=t,
                                                    tmax=tmax,
                                                    nrep=nrep,
                                                    c=c,
                                                    chopper_id=chopper_id,
                                                    bool_smooth=bool_smooth,
                                                    SG_w=SG_w,
                                                    SG_o=SG_o)
    id_ref = np.argmin(y0ref)

    id_i = np.shape(I)[0]
    id_j = np.shape(I)[1]
    id_t = np.shape(y0ref)[0]
    T_fobi = np.zeros((id_i, id_j, id_t))
    I0_fobi = np.zeros((id_i, id_j, id_t))
    I_fobi = np.zeros((id_i, id_j, id_t))
    for i in tqdm(range(0, id_i)):  #loop for all pixel position
        for j in range(0, id_j):
            [y0_rec, y_rec, T_rec,
             t_fobi] = full_fobi_reduction(I[i, j, :],
                                           I0[i, j, :],
                                           t=t,
                                           tmax=tmax,
                                           nrep=nrep,
                                           c=c,
                                           chopper_id=chopper_id,
                                           bool_smooth=bool_smooth,
                                           SG_w=SG_w,
                                           SG_o=SG_o)
            T_fobi[i, j, :] = np.roll(T_rec, -np.int(roll_value))
            I0_fobi[i, j, :] = np.roll(y0_rec, -np.int(roll_value))
            I_fobi[i, j, :] = np.roll(y_rec, -np.int(roll_value))

    return T_fobi, I0_fobi, I_fobi, t_fobi