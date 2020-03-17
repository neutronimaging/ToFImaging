import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.path as mplPath

def BOA_spectrum():
    BOA = np.loadtxt('C:\\Users\\busi_m\\Workspace\\neutronimaging\\ToFImaging\\lib\\dats\\BOA_spectrum.dat')
    return BOA
    
def ICON_spectrum():
    ICON = np.loadtxt('C:\\Users\\busi_m\\Workspace\\neutronimaging\\ToFImaging\\lib\\dats\\ICON_spectrum.dat')
    return ICON

def NEUTRA_spectrum():
    NEUTRA = np.loadtxt('C:\\Users\\busi_m\\Workspace\\neutronimaging\\ToFImaging\\lib\\dats\\NEUTRA_spectrum.dat')
    return NEUTRA

def POLDI_spectrum():
    POLDI = np.loadtxt('C:\\Users\\busi_m\\Workspace\\neutronimaging\\ToFImaging\\lib\\dats\\POLDI_spectrum.dat')
    return POLDI    