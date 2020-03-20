import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.path as mplPath

def BOA_spectrum():
    BOA = np.loadtxt('\\dats\\BOA_spectrum.dat')
    return BOA
    
def ICON_spectrum():
    ICON = np.loadtxt('\\ICON_spectrum.dat')
    return ICON

def NEUTRA_spectrum():
    NEUTRA = np.loadtxt('\\dats\\NEUTRA_spectrum.dat')
    return NEUTRA

def POLDI_spectrum():
    POLDI = np.loadtxt('\\dats\\POLDI_spectrum.dat')
    return POLDI    