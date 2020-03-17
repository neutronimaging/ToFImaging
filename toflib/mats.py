import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.path as mplPath

def NiO():
    NiO = np.loadtxt('C:\\Users\\busi_m\\Workspace\\neutronimaging\\ToFImaging\\lib\\dats\\NiO.dat')
    return NiO
    
def Fe_gamma():
    Fe_gamma = np.loadtxt('C:\\Users\\busi_m\\Workspace\\neutronimaging\\ToFImaging\\lib\\dats\\Fe_gamma.dat')
    return Fe_gamma

def Fe_alpha():
    Fe_alpha = np.loadtxt('C:\\Users\\busi_m\\Workspace\\neutronimaging\\ToFImaging\\lib\\dats\\Fe_alpha.dat')
    return Fe_alpha