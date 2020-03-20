import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.path as mplPath

def NiO():
    NiO = np.loadtxt('\\dats\\NiO.dat')
    return NiO
    
def Fe_gamma():
    Fe_gamma = np.loadtxt('dats\\Fe_gamma.dat')
    return Fe_gamma

def Fe_alpha():
    Fe_alpha = np.loadtxt('\\dats\\Fe_alpha.dat')
    return Fe_alpha