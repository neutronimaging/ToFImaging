import numpy as np
import TOF_routines

def chopper_time_delays_generator(time, nslits=8, nrep=2, mode='pseudorandom', rng = 0.25):
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

    if mode=='even':
        angles = np.linspace(0,360,nslits+1)/nrep
        angles = angles[1:-2]
        ang_t = angles
    if mode=='random':
        angles = 360*np.random.normal(size=nslits)/nrep
        ang_t = angles
    if mode=='pseudorandom':
        for i in range(0,nslits):
            angles[i] = i*360/(nslits*nrep) + np.random.normal()*360/(nslits*nrep)*rng
            ang_t = angles

    if nrep>1:
        for rep in range(0,nrep-1):
            angles = np.append(angles, ang_t+(rep+1)*360/nrep)
    
    shifts = Nt*angles/360
    D = np.zeros((Nt,1))
    for i in range(0,int(np.shape(shifts)[0])):
        sfloor = np.floor(shifts[i])
        rest = shifts[i]-sfloor
        sfloor = int(sfloor+1)
        D[sfloor] = 1-rest
        D[sfloor+1] = rest
    return D

def poldi_time_delays(time):
    """
    Generates time delays for the POLDI chopper:
    
    INPUTS:
    time = time-of-flight bins
    
    OUTPUT:
    D =  array with discretized dirac deltas at the chopper time delays.
    """
    Nt = int(np.shape(time)[0])
    angles = [0, 9.363, 21.475, 37.039, 50.417, 56.664, 67.422, 75.406]
    angles = [angles, angles+1*90, angles+2*90, angles+3*90]/360
    shifts = Nt*angles
    D = np.zeros((Nt,1))
    for i in range(0,np.shape(shifts)[0]):
        sfloor = np.floor(shifts[i])
        rest = shifts[i]-sfloor
        sfloor = sfloor+1
        D[sfloor] = 1-rest
        D[sfloor+1] = rest

    return D

def wiener_decorrelation(f, g, c=1e-2):
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
    arg = np.divide(F*G,(np.power(np.abs(G),2) + c))
    H = np.fft.ifft(arg)
    return H

def merge_reconstruction(x0,idx_l,nrep):
    """
    Merge repeated pattern signals into a single one
    
    INPUTS:
    x0 = repeated signal
    idx_l = index where the repeated signal starts
    nrep = number of pattern repetitions

    OUTPUT:
    x0_rec = merged signal
    """    
    nt = int(int(np.shape(x0)[0])/nrep)
    x0 = np.roll(x0,-idx_l)
    x0_rec = np.zeros(nt)
    for i in range(0,nrep,1):
        x0_rec = x0_rec + x0[i*nt:(i+1)*nt]
    x0_rec = x0_rec/nrep

    return x0_rec

def interp_noreadoutgaps(y,t,tmax,nrep=0,bool_plot=0):
    dt = np.nanmean(np.diff(t))
    t_tot = np.arange(t[0],t[-1]+dt,dt)
    app = np.nan*np.ones((np.shape(t_tot)[0]-np.shape(t)[0]))
    y[0] = np.nan
    y[-1] = np.nan
    y = np.concatenate((y,app))

    replen = np.floor(np.shape(y)[0]/nrep)
    y_overlap = np.zeros((replen,nrep))

    for i in range(0,nrep):
        y_overlap[:,0] = y[replen*(i-1):replen*i]

    y_merged = np.nanmean(y_overlap,axis=1)
    y_extended = y_merged
    for i in range(0,nrep-1):
        y_extended = np.concatenate((y_extended, y_merged))
    
    t_merged = t_tot[0:replen*nrep]

    return y_extended,t_merged

def full_fobi_reduction(y,y0,t,tmax,nrep,c=1e-1):
    [y,tn] = interp_noreadoutgaps(y,t,tmax,nrep)
    [y0,tn] = interp_noreadoutgaps(y0,t,tmax,nrep)

    D = poldi_time_delays(tn)

    x0rec = TOF_routines.savitzky_golay(wiener_decorrelation(y0,D,c))
    yrec = np.divide(TOF_routines.savitzky_golay(wiener_decorrelation(y,D,c)),x0rec)

    min_id = np.argmin(x0rec[0:np.shape(x0rec)[0]/nrep])
    x0rec = np.roll(x0rec,-min_id)
    yrec = np.roll(yrec,-min_id)

    replen = np.shape(y)[0]/nrep
    x0_over = np.zeros((replen,nrep))
    y_over = np.zeros((replen,nrep))
    for i in range(0,nrep):
        x0_over[:,i] = x0rec[replen*(i-1)+1:replen*i]
        y_over[:,i] = yrec[replen*(i-1)+1:replen*i]

    x0rec_merged = np.nanmean(x0_over,2)
    yrec_merged = np.nanmean(y_over,2)
    t_merged = tn[0:replen]

    return yrec_merged,x0rec_merged,t_merged