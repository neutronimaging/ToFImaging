import numpy as np
import tofimaging.ReductionTools as rt
from tqdm import tqdm


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


def interp_noreadoutgaps(y, t, tmax, nrep=0):
    """ Merge multislit chopper data interpolating readout gaps (in this case only at each chopper full rotation) then repeat the signal for the chopper pattern repetitions
    INPUTS:
    y = spectrum
    t = time-of-flight bins
    tmax = maximum time of flight (this parameter is dependent on the chopper frequency: tmax = 1/f  with f = chopper frequency)
    nrep = number of times the pattern is repeated (chopper)
    
    OUTPUTS:
    y_extended: interpolated signal
    t_merged: tof spectrum in the interpolated signal dimension
    """
    nrep = np.int(nrep)

    #dt = np.nanmean(np.diff(t),axis=0)
    #t_tot = np.arange(t[0],tmax+dt,dt)
    #app = np.nan*np.ones((np.shape(t_tot)[0]-np.shape(t)[0]))
    #y[0] = np.nan
    #y[-1] = np.nan
    #y_int = np.concatenate((y,app))
    #replen = np.int(np.floor(np.shape(y)[0]/nrep))

    replen = np.int(np.ceil(np.shape(y)[0] / nrep))
    t_tot = np.linspace(t[0], tmax, nrep * replen)
    y_int = np.interp(t_tot, t, y, left=np.nan, right=np.nan)

    y_overlap = np.zeros((replen, nrep))
    for i in range(0, nrep):
        y_overlap[:, i] = y_int[replen * i:replen * (i + 1)]

    y_merged = np.nanmedian(y_overlap, axis=1)
    y_extended = y_merged
    for i in range(0, nrep - 1):
        y_extended = np.concatenate((y_extended, y_merged))

    t_extended = t_tot[0:replen * nrep]

    return y_extended, t_extended


def full_fobi_reduction(y,
                        y0,
                        t,
                        tmax,
                        nrep,
                        chopper_id,
                        c=1e-1,
                        bool_roll=False,
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
    [y, tn] = interp_noreadoutgaps(y, t, tmax, nrep)
    [y0, tn] = interp_noreadoutgaps(y0, t, tmax, nrep)
    if (chopper_id == 'poldi'):
        D = time_delays_poldi(tn)
    if (chopper_id == '4x10'):
        D = time_delays_4x10(tn)
    if (chopper_id == '5x8'):
        D = time_delays_5x8(tn)

    # FIXED POLDI angles: now it's working as deconvolution!
    # x0rec = rt.savitzky_golay(wiener_decorrelation(y0,D,c),SG_w,SG_o)
    # yrec = rt.savitzky_golay(wiener_decorrelation(y,D,c),SG_w,SG_o)
    # Trec = np.divide(yrec,x0rec)
    # #Trec = rt.savitzky_golay(wiener_decorrelation(np.divide(y,y0),D,c),SG_w,SG_o)

    y0rec = wiener_deconvolution(y0, D, c)
    yrec = wiener_deconvolution(y, D, c)
    if (bool_smooth):
        y0rec = rt.savitzky_golay(y0rec, SG_w, SG_o)
        yrec = rt.savitzky_golay(yrec, SG_w, SG_o)
    Trec = np.divide(yrec, y0rec)
    # Trec = rt.savitzky_golay(wiener_deconvolution(np.divide(y,y0),D,c),SG_w,SG_o)

    if (bool_roll == True):
        min_id = np.argmin(y0rec[0:np.shape(y0rec)[0] / nrep])
        y0rec = np.roll(y0rec, -min_id)
        yrec = np.roll(yrec, -min_id)
        Trec = np.roll(Trec, -min_id)

    replen = np.int(np.floor(np.shape(y)[0] / nrep))
    y0_over = np.zeros((replen, nrep))
    y_over = np.zeros((replen, nrep))
    T_over = np.zeros((replen, nrep))
    for i in range(0, nrep):
        y0_over[:, i] = y0rec[replen * i:replen * (i + 1)]
        y_over[:, i] = yrec[replen * i:replen * (i + 1)]
        T_over[:, i] = Trec[replen * i:replen * (i + 1)]

    y0_fobi = np.nanmean(y0_over, axis=1)
    y_fobi = np.nanmean(y_over, axis=1)
    T_fobi = np.nanmean(T_over, axis=1)
    t_fobi = tn[0:replen]

    return y0_fobi, y_fobi, T_fobi, t_fobi


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