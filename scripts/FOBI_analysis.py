import numpy as np

def chopper_time_delays(time, nslits=8, nrep=2, mode='pseudorandom', rng = 0.25):
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

def wiener_decorrelation(f, g, c=1e-2):
    F = np.fft.fft(f) 
    G = np.fft.fft(g)
    arg = np.divide(F*G,(np.power(np.abs(G),2) + c))
    H = np.fft.ifft(arg)
    return H

def merge_reconstruction(x0,idx_l,nrep):
    nt = int(int(np.shape(x0)[0])/nrep)
    x0 = np.roll(x0,-idx_l)
    x0_rec = np.zeros(nt)
    for i in range(0,nrep,1):
        x0_rec = x0_rec + x0[i*nt:(i+1)*nt]
    x0_rec = x0_rec/nrep

    return x0_rec