import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.path as mplPath

def rmse(y,y0):
    T = int(np.shape(y)[0])
    rmse = np.sqrt(np.sum(np.power(y-y0,2))/T)
    return rmse