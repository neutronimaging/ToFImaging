import numpy as np
from numpy import pi, r_, math, random
import matplotlib.pyplot as plt
from scipy import special
from lmfit import Model
import math

from reduction_tools import find_nearest
from reduction_tools import savitzky_golay as SG_filter

def phase_ratio_linearcomb(lac,spectrum,phase1lac,phase2lac,phase_spectrum,lambda_range_norm,lambda_range_fit,est_phi=0.5,method ='least_squares',bool_SG=False,SG_w=5,SG_n=1,bool_print=False): 
    """ Performs phase ratio fitting on linear combination of two basis functions, works with linear attenuation coefficient (LAC) spectra
    INPUTS:
    lac = 1d array the attenuation -log(I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
    phase1lac = lac of the phase 1
    phase2lac = lac of the phase 2
    spectrum_phase = spectrum corresponding to phase 1 and 2
    lambda_range_norm = lambda range where to normalize spectra
    lambda_range_edges = lambda range where to do the fitting
    est_phi = estimate phase 1 weight
    method = fitting method
    bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative
    SG_w = window size of S-G filter
    SG_n = order of S-G filter
    bool_print = set to True to print output

    OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    'phi' : phase1 weight
    """ 

    if(bool_SG):
        lac = SG_filter(lac,SG_w,SG_n)
    lac0 = lac
    sp0 = spectrum
    #cutting to selected range for normalization
    idx_low = find_nearest(spectrum,lambda_range_norm[0])
    idx_high = find_nearest(spectrum,lambda_range_norm[1])
    lac = lac[idx_low:idx_high]
    spectrum = spectrum[idx_low:idx_high]

    idx_low = find_nearest(phase_spectrum,lambda_range_norm[0])
    idx_high = find_nearest(phase_spectrum,lambda_range_norm[1])
    lambda_table = phase_spectrum[idx_low:idx_high]
    phase1lac = phase1lac[idx_low:idx_high]
    phase2lac = phase2lac[idx_low:idx_high]

    phase1lac = np.interp(spectrum,lambda_table,phase1lac)
    phase2lac = np.interp(spectrum,lambda_table,phase2lac)
    phase1lac = phase1lac*(np.nanmean(lac)/np.nanmean(phase1lac))
    phase2lac = phase2lac*(np.nanmean(lac)/np.nanmean(phase2lac))

    # initial conditions
    idx_l = find_nearest(spectrum,lambda_range_fit[0])
    idx_h = find_nearest(spectrum,lambda_range_fit[1])
    lac_cut = lac[idx_l:idx_h]
    ph1 = phase1lac[idx_l:idx_h]
    ph2 = phase2lac[idx_l:idx_h]
    l_cut = spectrum[idx_l:idx_h]

    # method='trust_exact'
    # method='nelder' #not bad
    # method='differential_evolution' # needs bounds
    # method='basinhopping' # not bad
    # method='lmsquare' # this should implement the Levemberq-Marquardt but it says Nelder-Mead method (which should be Amoeba)
    # method ='least_squares' # default and it implements the Levenberg-Marquardt

    def two_phases(ph1,ph2,f):
        return f*ph1+(1-f)*ph2

    gmodel = Model(two_phases,independent_vars=['ph1', 'ph2'])
    params = gmodel.make_params(f = est_phi)    
    #params['ph1'].vary = False
    #params['ph2'].vary = False
    params['f'].vary = True
    params['f'].min = 0.0
    params['f'].max = 1.0
    result = gmodel.fit(lac_cut, params, ph1 = ph1, ph2 = ph2, method=method, nan_policy='propagate')
    phi=result.best_values.get('f')    

    if(bool_print):
        print('phase fraction (ph1 %) = ',100*phi,'%')
        plt.figure()
        plt.plot(sp0, lac0, label ='LAC')
        plt.plot(spectrum, lac, label ='LAC (normalization)')
        plt.plot(spectrum, phase1lac, '--', label ='Phase 1')
        plt.plot(spectrum, phase2lac, '--', label ='Phase 2')
        plt.plot(l_cut, two_phases(ph1,ph2,phi), label ='Fitted ratio')
        plt.title('Bragg pattern'), plt.xlabel('Wavelenght [Å]')
        plt.legend(),        plt.show(),        plt.close()

    return {'phi': phi}

def phase_ratio_linearcomb_three(lac,spectrum,phase1lac,phase2lac,phase3lac,phase_spectrum,lambda_range_norm,lambda_range_fit,est_f1=0.333,est_f2=0.333,est_f3=0.334,method ='least_squares',bool_SG=False,SG_w=5,SG_n=1,bool_print=False):     
    """ Performs phase ratio fitting on linear combination of two basis functions, works with linear attenuation coefficient (LAC) spectra
    INPUTS:
    lac = 1d array the attenuation -log(I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
    phase1lac = lac of the phase 1
    phase2lac = lac of the phase 2
    phase3lac = lac of the phase 3
    spectrum_phase = spectrum corresponding to phase 1 and 2
    lambda_range_norm = lambda range where to normalize spectra
    lambda_range_edges = lambda range where to do the fitting
    est_f1 = estimate phase 1 weight
    est_f2 = estimate phase 2 weight
    est_f3 = estimate phase 3 weight
    method = fitting method
    bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative
    SG_w = window size of S-G filter
    SG_n = order of S-G filter
    bool_print = set to True to print output

    OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    'phi1' : phase1 weight
    'phi2' : phase2 weight
    'phi3' : phase3 weight
    """ 
    
    if(bool_SG):
        lac = SG_filter(lac,SG_w,SG_n)
    lac0 = lac
    sp0 = spectrum
    #cutting to selected range for normalization
    idx_low = find_nearest(spectrum,lambda_range_norm[0])
    idx_high = find_nearest(spectrum,lambda_range_norm[1])
    lac = lac[idx_low:idx_high]
    spectrum = spectrum[idx_low:idx_high]

    idx_low = find_nearest(phase_spectrum,lambda_range_norm[0])
    idx_high = find_nearest(phase_spectrum,lambda_range_norm[1])
    lambda_table = phase_spectrum[idx_low:idx_high]
    phase1lac = phase1lac[idx_low:idx_high]
    phase2lac = phase2lac[idx_low:idx_high]
    phase3lac = phase3lac[idx_low:idx_high]

    phase1lac = np.interp(spectrum,lambda_table,phase1lac)
    phase2lac =  np.interp(spectrum,lambda_table,phase2lac)
    phase3lac =  np.interp(spectrum,lambda_table,phase3lac)
    phase1lac = phase1lac*(np.nanmean(lac)/np.nanmean(phase1lac))
    phase2lac = phase2lac*(np.nanmean(lac)/np.nanmean(phase2lac))
    phase3lac = phase3lac*(np.nanmean(lac)/np.nanmean(phase3lac))

    # initial conditions
    idx_l = find_nearest(spectrum,lambda_range_fit[0])
    idx_h = find_nearest(spectrum,lambda_range_fit[1])
    lac_cut = lac[idx_l:idx_h]
    ph1 = phase1lac[idx_l:idx_h]
    ph2 = phase2lac[idx_l:idx_h]
    ph3 = phase3lac[idx_l:idx_h]
    l_cut = spectrum[idx_l:idx_h]

    # method='trust_exact'
    # method='nelder' #not bad
    # method='differential_evolution' # needs bounds
    # method='basinhopping' # not bad
    # method='lmsquare' # this should implement the Levemberq-Marquardt but it says Nelder-Mead method (which should be Amoeba)
    # method ='least_squares' # default and it implements the Levenberg-Marquardt

    def three_phases(ph1,ph2,ph3,f1,f2,f3):
        return f1*ph1+f2*ph2+f3*ph3
        # return (1-f1-f2)*ph1+f1*ph2+f2*ph3

    gmodel = Model(three_phases,independent_vars=['ph1', 'ph2', 'ph3'])
    params = gmodel.make_params(f1 = est_f1, f2 = est_f2, f3 = est_f3)    
    #params['ph1'].vary = False
    #params['ph2'].vary = False
    params['f1'].vary = True
    params['f1'].min = 0.0
    params['f1'].max = 1.0
    params['f2'].vary = True
    params['f2'].min = 0.0
    params['f2'].max = 1.0
    params['f3'].vary = True
    params['f3'].min = 0.0
    params['f3'].max = 1.0
    result = gmodel.fit(lac_cut, params, ph1 = ph1, ph2 = ph2, ph3 = ph3, method=method, nan_policy='propagate')
    phi1=result.best_values.get('f1')    
    phi2=result.best_values.get('f2')    
    phi3=result.best_values.get('f3')    
    phi1_n = phi1/(phi1+phi2+phi3)
    phi2_n = phi2/(phi1+phi2+phi3)
    phi3_n = phi3/(phi1+phi2+phi3)

    if(bool_print):
        print('Phase fraction 1 (ph1 %) = ',100*phi1,'%')
        print('Phase fraction 2 (ph2 %) = ',100*phi2,'%')
        print('Phase fraction 3 (ph3 %) = ',100*phi3,'%')
        plt.figure()
        plt.plot(sp0, lac0, label ='LAC')
        plt.plot(spectrum, lac, label ='LAC (normalization)')
        plt.plot(spectrum, phase1lac, '--', label ='Phase 1')
        plt.plot(spectrum, phase2lac, '--', label ='Phase 2')
        plt.plot(spectrum, phase3lac, '--', label ='Phase 3')
        plt.plot(l_cut, three_phases(ph1,ph2,ph3,phi1,phi2,phi3), label ='Fitted ratio')
        plt.title('Bragg pattern'), plt.xlabel('Wavelenght [Å]')
        plt.legend(),        plt.show(),        plt.close()

    return {'phi1': phi1_n,'phi2': phi2_n,'phi3': phi3_n}