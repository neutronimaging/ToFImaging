import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from lmfit import Model
from tqdm import tqdm

import tofimaging.ReductionTools as rt
from tofimaging.ReductionTools import find_nearest
from tofimaging.ReductionTools import savitzky_golay as SG_filter


def PhaseRatioLinearCombination(lac,
                                spectrum,
                                phase1lac,
                                phase2lac,
                                phase_spectrum,
                                lambda_range_norm,
                                lambda_range_fit,
                                est_phi=0.5,
                                method='least_squares',
                                bool_SG=False,
                                SG_w=5,
                                SG_n=1,
                                bool_print=False):
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

    if (bool_SG):
        lac = SG_filter(lac, SG_w, SG_n)
    lac0 = lac
    sp0 = spectrum
    #cutting to selected range for normalization
    idx_low = find_nearest(spectrum, lambda_range_norm[0])
    idx_high = find_nearest(spectrum, lambda_range_norm[1])
    lac = lac[idx_low:idx_high]
    spectrum = spectrum[idx_low:idx_high]

    idx_low = find_nearest(phase_spectrum, lambda_range_norm[0])
    idx_high = find_nearest(phase_spectrum, lambda_range_norm[1])
    lambda_table = phase_spectrum[idx_low:idx_high]
    phase1lac = phase1lac[idx_low:idx_high]
    phase2lac = phase2lac[idx_low:idx_high]

    phase1lac = np.interp(spectrum, lambda_table, phase1lac)
    phase2lac = np.interp(spectrum, lambda_table, phase2lac)
    phase1lac = phase1lac * (np.nanmean(lac) / np.nanmean(phase1lac))
    phase2lac = phase2lac * (np.nanmean(lac) / np.nanmean(phase2lac))

    # initial conditions
    idx_l = find_nearest(spectrum, lambda_range_fit[0])
    idx_h = find_nearest(spectrum, lambda_range_fit[1])
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

    def two_phases(ph1, ph2, f):
        return f * ph1 + (1 - f) * ph2

    gmodel = Model(two_phases, independent_vars=['ph1', 'ph2'])
    params = gmodel.make_params(f=est_phi)
    #params['ph1'].vary = False
    #params['ph2'].vary = False
    params['f'].vary = True
    params['f'].min = 0.0
    params['f'].max = 1.0
    result = gmodel.fit(lac_cut,
                        params,
                        ph1=ph1,
                        ph2=ph2,
                        method=method,
                        nan_policy='propagate')
    phi = result.best_values.get('f')

    if (bool_print):
        print('phase fraction (ph1 %) = ', 100 * phi, '%')
        plt.figure()
        plt.plot(sp0, lac0, label='LAC')
        plt.plot(spectrum, lac, label='LAC (normalization)')
        plt.plot(spectrum, phase1lac, '--', label='Phase 1')
        plt.plot(spectrum, phase2lac, '--', label='Phase 2')
        plt.plot(l_cut, two_phases(ph1, ph2, phi), label='Fitted ratio')
        plt.title('Bragg pattern'), plt.xlabel('Wavelenght [Å]')
        plt.legend(), plt.show(), plt.close()

    return {'phi': phi}


def PhaseRatioLinearCombination2D(lac_tof,
                                  spectrum,
                                  phase1lac,
                                  phase2lac,
                                  spectrum_phase,
                                  lambda_range_norm,
                                  lambda_range_edges,
                                  calibration_matrix=np.ndarray([0]),
                                  mask=np.ndarray([0]),
                                  auto_mask=True,
                                  mask_thresh=[0.05, 0.95],
                                  est_phi=0.5,
                                  method='least_squares',
                                  bool_SG=False,
                                  SG_w=5,
                                  SG_n=1,
                                  bool_save=False,
                                  bool_print=False,
                                  debug_idx=[]):
    """ Performs phase ratio fitting on linear combination of two basis functions, works with linear attenuation coefficient (LAC) spectra
    INPUTS:
    lac_tof = 3d matrix with the stack of attenuation -log(I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
    phase1lac = lac of the phase 1
    phase2lac = lac of the phase 2
    spectrum_phase = spectrum corresponding to phase 1 and 2
    lambda_range_norm = lambda range where to normalize spectra
    lambda_range_edges = lambda range where to do the fitting
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y)
    mask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_phi = estimated phase fraction
    method = fitting method
    bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative
    SG_w = window size of S-G filter
    SG_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    dictionary with the following fit in the dimension of the mask
    'phase_ratio' : phase ratio
    """

    if (mask.any()):
        mymask = mask
        plt.figure()
        plt.subplot(1, 2, 1), plt.imshow(np.median(
            lac_tof, axis=2)), plt.title('Full-spectrum Image')
        plt.subplot(1, 2, 2), plt.imshow(mymask), plt.title('Mask')
        plt.tight_layout()
        plt.show()
        plt.close()
        if ([np.shape(lac_tof)[0], np.shape(lac_tof)[1]] !=
            [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif (auto_mask):
        import skimage.filters
        mymask = rt.medianimage(lac_tof)
        plt.figure()
        plt.subplot(1, 3,
                    1), plt.imshow(mymask), plt.title('Full-spectrum Image')
        mymask[mymask > mask_thresh[1]] = 0.0
        mymask[mymask < mask_thresh[0]] = 0.0
        mymask[mymask > 0] = 1.0
        mymask[np.isinf(mymask)] = 0.0
        mymask[np.isnan(mymask)] = 0.0
        plt.subplot(1, 3, 2), plt.imshow(mymask), plt.title('Mask')
        mymask = skimage.filters.gaussian(mymask, sigma=2)
        mymask[mymask > 0] = 1.0
        plt.subplot(1, 3, 3), plt.imshow(mymask), plt.title('Mask - gauss')
        plt.tight_layout()
        plt.show()
        plt.close()
    else:
        mymask = np.ones([np.shape(lac_tof)[0], np.shape(lac_tof)[1]])

    if (calibration_matrix.any()):
        if ((np.shape(lac_tof)[0] != np.shape(calibration_matrix)[0]) |
            (np.shape(lac_tof)[1] != np.shape(calibration_matrix)[1])):
            print(
                '!!!!!!! WARNING CALIBRATION MATRIX HAS NOT SAME SIZE OF IMAGE !!!!!!!!!!!!!!'
            )

    if (any(debug_idx)):  #testing on a single pixel
        lac = lac_tof[debug_idx[0], debug_idx[1], :]

        if (calibration_matrix.any()):
            lambd = rt.tof2l_t0k(
                spectrum, calibration_matrix[debug_idx[0], debug_idx[1], 1],
                calibration_matrix[debug_idx[0], debug_idx[1], 0])
        else:
            lambd = spectrum

        PhaseRatioLinearCombination(lac,
                                    lambd,
                                    phase1lac,
                                    phase2lac,
                                    spectrum_phase,
                                    lambda_range_norm,
                                    lambda_range_edges,
                                    est_phi=est_phi,
                                    method=method,
                                    bool_SG=bool_SG,
                                    SG_w=SG_w,
                                    SG_n=SG_n,
                                    bool_print=1)
        return

    phase_ratio = np.zeros(np.shape(mymask))
    for i in tqdm(range(0, np.shape(mymask)[0])):
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i, j]):
                lac = lac_tof[i, j, :]

                if (calibration_matrix.any()):
                    lambd = rt.tof2l_t0k(spectrum, calibration_matrix[i, j, 1],
                                         calibration_matrix[i, j, 0])
                else:
                    lambd = spectrum

                try:
                    phi_fit = PhaseRatioLinearCombination(lac,
                                                          lambd,
                                                          phase1lac,
                                                          phase2lac,
                                                          spectrum_phase,
                                                          lambda_range_norm,
                                                          lambda_range_edges,
                                                          est_phi=est_phi,
                                                          method=method,
                                                          bool_SG=bool_SG,
                                                          SG_w=SG_w,
                                                          SG_n=SG_n,
                                                          bool_print=0)
                    phase_ratio[i, j] = phi_fit['phi']
                except:
                    print("Unexpected error at :", i, j)
                    phase_ratio[i, j] = -2.0

    if (bool_print):
        plt.figure()
        plt.imshow(phase_ratio, cmap='jet'), plt.title('Phase ratio (%)')
        plt.colorbar()
        plt.show(), plt.close()
    if (bool_save):
        np.save('phase_ratio.npy', phase_ratio)
    return {'phase_ratio': phase_ratio}


def PhaseRatioLinearCombination3(lac,
                                 spectrum,
                                 phase1lac,
                                 phase2lac,
                                 phase3lac,
                                 phase_spectrum,
                                 lambda_range_norm,
                                 lambda_range_fit,
                                 est_f1=0.333,
                                 est_f2=0.333,
                                 est_f3=0.334,
                                 method='least_squares',
                                 bool_SG=False,
                                 SG_w=5,
                                 SG_n=1,
                                 bool_print=False):
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

    if (bool_SG):
        lac = SG_filter(lac, SG_w, SG_n)
    lac0 = lac
    sp0 = spectrum
    #cutting to selected range for normalization
    idx_low = find_nearest(spectrum, lambda_range_norm[0])
    idx_high = find_nearest(spectrum, lambda_range_norm[1])
    lac = lac[idx_low:idx_high]
    spectrum = spectrum[idx_low:idx_high]

    idx_low = find_nearest(phase_spectrum, lambda_range_norm[0])
    idx_high = find_nearest(phase_spectrum, lambda_range_norm[1])
    lambda_table = phase_spectrum[idx_low:idx_high]
    phase1lac = phase1lac[idx_low:idx_high]
    phase2lac = phase2lac[idx_low:idx_high]
    phase3lac = phase3lac[idx_low:idx_high]

    phase1lac = np.interp(spectrum, lambda_table, phase1lac)
    phase2lac = np.interp(spectrum, lambda_table, phase2lac)
    phase3lac = np.interp(spectrum, lambda_table, phase3lac)
    phase1lac = phase1lac * (np.nanmean(lac) / np.nanmean(phase1lac))
    phase2lac = phase2lac * (np.nanmean(lac) / np.nanmean(phase2lac))
    phase3lac = phase3lac * (np.nanmean(lac) / np.nanmean(phase3lac))

    # initial conditions
    idx_l = find_nearest(spectrum, lambda_range_fit[0])
    idx_h = find_nearest(spectrum, lambda_range_fit[1])
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

    def three_phases(ph1, ph2, ph3, f1, f2, f3):
        return f1 * ph1 + f2 * ph2 + f3 * ph3
        # return (1-f1-f2)*ph1+f1*ph2+f2*ph3

    gmodel = Model(three_phases, independent_vars=['ph1', 'ph2', 'ph3'])
    params = gmodel.make_params(f1=est_f1, f2=est_f2, f3=est_f3)
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
    result = gmodel.fit(lac_cut,
                        params,
                        ph1=ph1,
                        ph2=ph2,
                        ph3=ph3,
                        method=method,
                        nan_policy='propagate')
    phi1 = result.best_values.get('f1')
    phi2 = result.best_values.get('f2')
    phi3 = result.best_values.get('f3')
    phi1_n = phi1 / (phi1 + phi2 + phi3)
    phi2_n = phi2 / (phi1 + phi2 + phi3)
    phi3_n = phi3 / (phi1 + phi2 + phi3)

    if (bool_print):
        print('Phase fraction 1 (ph1 %) = ', 100 * phi1, '%')
        print('Phase fraction 2 (ph2 %) = ', 100 * phi2, '%')
        print('Phase fraction 3 (ph3 %) = ', 100 * phi3, '%')
        plt.figure()
        plt.plot(sp0, lac0, label='LAC')
        plt.plot(spectrum, lac, label='LAC (normalization)')
        plt.plot(spectrum, phase1lac, '--', label='Phase 1')
        plt.plot(spectrum, phase2lac, '--', label='Phase 2')
        plt.plot(spectrum, phase3lac, '--', label='Phase 3')
        plt.plot(l_cut,
                 three_phases(ph1, ph2, ph3, phi1, phi2, phi3),
                 label='Fitted ratio')
        plt.title('Bragg pattern'), plt.xlabel('Wavelenght [Å]')
        plt.legend(), plt.show(), plt.close()

    return {'phi1': phi1_n, 'phi2': phi2_n, 'phi3': phi3_n}


def PhaseRatioLinearCombination32D(lac_tof,
                                   spectrum,
                                   phase1lac,
                                   phase2lac,
                                   phase3lac,
                                   spectrum_phase,
                                   lambda_range_norm,
                                   lambda_range_edges,
                                   calibration_matrix=np.ndarray([0]),
                                   mask=np.ndarray([0]),
                                   auto_mask=True,
                                   mask_thresh=[0.05, 0.95],
                                   est_f1=0.333,
                                   est_f2=0.333,
                                   est_f3=0.334,
                                   method='least_squares',
                                   bool_SG=False,
                                   SG_w=5,
                                   SG_n=1,
                                   bool_save=False,
                                   bool_print=False,
                                   debug_idx=[]):
    """ Performs phase ratio fitting on linear combination of two basis functions, works with linear attenuation coefficient (LAC) spectra
    INPUTS:
    lac_tof = 3d matrix with the stack of attenuation -log(I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
    phase1lac = lac of the phase 1
    phase2lac = lac of the phase 2
    spectrum_phase = spectrum corresponding to phase 1 and 2
    lambda_range_norm = lambda range where to normalize spectra
    lambda_range_edges = lambda range where to do the fitting
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y)
    mask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    est_f1 = estimate phase 1 weight
    est_f2 = estimate phase 2 weight
    est_f3 = estimate phase 3 weight
    method = fitting method
    bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative
    SG_w = window size of S-G filter
    SG_n = order of S-G filter
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    dictionary with the following fit in the dimension of the mask
    'phase1_ratio' : phase1 weight
    'phase2_ratio' : phase2 weight
    'phase3_ratio' : phase3 weight
    """

    if (mask.any()):
        mymask = mask
        plt.figure()
        plt.subplot(1, 2, 1), plt.imshow(np.median(
            lac_tof, axis=2)), plt.title('Full-spectrum Image')
        plt.subplot(1, 2, 2), plt.imshow(mymask)
        plt.title('Mask')
        plt.tight_layout()
        plt.show()
        plt.close()
        if ([np.shape(lac_tof)[0], np.shape(lac_tof)[1]] !=
            [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif (auto_mask):
        import skimage.filters
        mymask = rt.medianimage(lac_tof)
        plt.figure()
        plt.subplot(1, 3,
                    1), plt.imshow(mymask), plt.title('Full-spectrum Image')
        mymask[mymask > mask_thresh[1]] = 0.0
        mymask[mymask < mask_thresh[0]] = 0.0
        mymask[mymask > 0] = 1.0
        mymask[np.isinf(mymask)] = 0.0
        mymask[np.isnan(mymask)] = 0.0
        plt.subplot(1, 3, 2), plt.imshow(mymask), plt.title('Mask')
        mymask = skimage.filters.gaussian(mymask, sigma=2)
        mymask[mymask > 0] = 1.0
        plt.subplot(1, 3, 3), plt.imshow(mymask), plt.title('Mask - gauss')
        plt.tight_layout()
        plt.show()
        plt.close()
    else:
        mymask = np.ones([np.shape(lac_tof)[0], np.shape(lac_tof)[1]])

    if (calibration_matrix.any()):
        if ((np.shape(lac_tof)[0] != np.shape(calibration_matrix)[0]) |
            (np.shape(lac_tof)[1] != np.shape(calibration_matrix)[1])):
            print(
                '!!!!!!! WARNING CALIBRATION MATRIX HAS NOT SAME SIZE OF IMAGE !!!!!!!!!!!!!!'
            )

    if (any(debug_idx)):  #testing on a single pixel
        lac = lac_tof[debug_idx[0], debug_idx[1], :]

        if (calibration_matrix.any()):
            lambd = rt.tof2l_t0k(
                spectrum, calibration_matrix[debug_idx[0], debug_idx[1], 1],
                calibration_matrix[debug_idx[0], debug_idx[1], 0])
        else:
            lambd = spectrum

        PhaseRatioLinearCombination3(lac,
                                     lambd,
                                     phase1lac,
                                     phase2lac,
                                     phase3lac,
                                     spectrum_phase,
                                     lambda_range_norm,
                                     lambda_range_edges,
                                     est_f1=est_f1,
                                     est_f2=est_f2,
                                     method=method,
                                     bool_SG=bool_SG,
                                     SG_w=SG_w,
                                     SG_n=SG_n,
                                     bool_print=1)
        return

    phase1_ratio = np.zeros(np.shape(mymask))
    phase2_ratio = np.zeros(np.shape(mymask))
    phase3_ratio = np.zeros(np.shape(mymask))
    for i in tqdm(range(0, np.shape(mymask)[0])):
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i, j]):
                lac = lac_tof[i, j, :]

                if (calibration_matrix.any()):
                    lambd = rt.tof2l_t0k(spectrum, calibration_matrix[i, j, 1],
                                         calibration_matrix[i, j, 0])
                else:
                    lambd = spectrum

                try:
                    phi_fit = PhaseRatioLinearCombination3(lac,
                                                           lambd,
                                                           phase1lac,
                                                           phase2lac,
                                                           phase3lac,
                                                           spectrum_phase,
                                                           lambda_range_norm,
                                                           lambda_range_edges,
                                                           est_f1=est_f1,
                                                           est_f2=est_f2,
                                                           method=method,
                                                           bool_SG=bool_SG,
                                                           SG_w=SG_w,
                                                           SG_n=SG_n,
                                                           bool_print=0)
                    phase1_ratio[i, j] = phi_fit['phi1']
                    phase2_ratio[i, j] = phi_fit['phi2']
                    phase3_ratio[i, j] = phi_fit['phi3']
                except:
                    print("Unexpected error at :", i, j)
                    phase1_ratio[i, j] = -2.0
                    phase2_ratio[i, j] = -2.0
                    phase3_ratio[i, j] = -2.0

    if (bool_print):
        plt.figure()
        plt.subplot(1, 3, 1), plt.imshow(
            phase1_ratio,
            cmap='jet'), plt.title('Phase 1 weight (%)'), plt.colorbar()
        plt.subplot(1, 3, 2), plt.imshow(
            phase2_ratio,
            cmap='jet'), plt.title('Phase 2 weight (%)'), plt.colorbar()
        plt.subplot(1, 3, 3), plt.imshow(
            phase3_ratio,
            cmap='jet'), plt.title('Phase 3 weight (%)'), plt.colorbar()
        plt.tight_layout()
        plt.show()
        plt.close()
    if (bool_save):
        np.save('phase1_ratio.npy', phase1_ratio)
        np.save('phase2_ratio.npy', phase2_ratio)
        np.save('phase3_ratio.npy', phase3_ratio)
    return {
        'phase1_ratio': phase1_ratio,
        'phase2_ratio': phase2_ratio,
        'phase3_ratio': phase3_ratio
    }


def WavelengthSelectiveRatio(lac,
                             spectrum,
                             l1,
                             l2,
                             l1_w=0,
                             l2_w=0,
                             bool_SG=False,
                             SG_w=5,
                             SG_n=1,
                             bool_print=False):
    """ Performs phase ratio fitting on linear combination of two basis functions, works with linear attenuation coefficient (LAC) spectra
    INPUTS:
    lac = 1d array the attenuation -log(I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
    l1 = lambda position of the 1st point
    l2 = lambda position of the 2nd point
    l1_w = lambda window size (bi-directional) of the 1st point
    l2_w = lambda window size (bi-directional) of the 2nd point
    bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative
    SG_w = window size of S-G filter
    SG_n = order of S-G filter

    OUTPUTS:
    #dictionary with the following fit in the dimension of the mask
    'WSR' : Wavelength Selective Ratio
    """

    if (bool_SG):
        lac = SG_filter(lac, SG_w, SG_n)

    idx_1 = find_nearest(spectrum, l1)
    idx_2 = find_nearest(spectrum, l2)

    WSR = np.nanmean(lac[idx_1 - l1_w:idx_1 + l1_w + 1]) / np.nanmean(
        lac[idx_2 - l2_w:idx_2 + l2_w + 1])

    if (bool_print):
        plt.plot(spectrum, lac)
        plt.plot(spectrum[idx_1 - l1_w:idx_1 + l1_w + 1],
                 lac[idx_1 - l1_w:idx_1 + l1_w + 1],
                 'x',
                 markeredgewidth=3,
                 c='orange')
        plt.plot(spectrum[idx_2 - l2_w:idx_2 + l2_w + 1],
                 lac[idx_2 - l2_w:idx_2 + l2_w + 1],
                 'o',
                 markeredgewidth=3,
                 c='orange')
        print(WSR)

    return {'WSR': WSR}


def WavelengthSelectiveRatio2D(lac_tof,
                               spectrum,
                               l1,
                               l2,
                               l1_w=0,
                               l2_w=0,
                               calibration_matrix=np.ndarray([0]),
                               mask=np.ndarray([0]),
                               auto_mask=True,
                               mask_thresh=[0.05, 0.95],
                               bool_SG=False,
                               SG_w=5,
                               SG_n=1,
                               bool_save=False,
                               bool_print=False,
                               debug_idx=[]):
    """ Performs phase ratio fitting on linear combination of two basis functions, works with linear attenuation coefficient (LAC) spectra
    INPUTS:
    lac_tof = 3d matrix with the stack of attenuation -log(I/I0) TOF images (x,y,lambda) 
    spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
    l1 = lambda position of the 1st point
    l2 = lambda position of the 2nd point
    l1_w = lambda window size (bi-directional) of the 1st point
    l2_w = lambda window size (bi-directional) of the 2nd point
    calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]);
                         will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y)
    mask = mask of where to perform the fit (x,y)
    auto_mask = if True will automatically mask the region based on the mask_thresh thresholds
    mask_thresh = low and high threshold for the automatic mask
    bool_save = set to True to save output
    bool_print = set to True to print output
    debug_idx = pixel coordinates where to test the single pixel fitting

    OUTPUTS:
    dictionary with the following fit in the dimension of the mask
    'WSR' : Wavelength Selective Ratio
    """

    if (mask.any()):
        mymask = mask
        plt.figure()
        plt.subplot(1, 2, 1), plt.imshow(np.median(
            lac_tof, axis=2)), plt.title('Full-spectrum Image')
        plt.subplot(1, 2, 2), plt.imshow(mymask), plt.title('Mask')
        plt.tight_layout()
        plt.show()
        plt.close()
        if ([np.shape(lac_tof)[0], np.shape(lac_tof)[1]] !=
            [np.shape(mymask)[0], np.shape(mymask)[1]]):
            print('WARNING: Mask size does not match frames size')
    elif (auto_mask):
        import skimage.filters
        mymask = rt.medianimage(lac_tof)
        plt.figure()
        plt.subplot(1, 3,
                    1), plt.imshow(mymask), plt.title('Full-spectrum Image')
        mymask[mymask > mask_thresh[1]] = 0.0
        mymask[mymask < mask_thresh[0]] = 0.0
        mymask[mymask > 0] = 1.0
        mymask[np.isinf(mymask)] = 0.0
        mymask[np.isnan(mymask)] = 0.0
        plt.subplot(1, 3, 2), plt.imshow(mymask), plt.title('Mask')
        mymask = skimage.filters.gaussian(mymask, sigma=2)
        mymask[mymask > 0] = 1.0
        plt.subplot(1, 3, 3), plt.imshow(mymask), plt.title('Mask - gauss')
        plt.tight_layout()
        plt.show()
        plt.close()
    else:
        mymask = np.ones([np.shape(lac_tof)[0], np.shape(lac_tof)[1]])

    if (calibration_matrix.any()):
        if ((np.shape(lac_tof)[0] != np.shape(calibration_matrix)[0]) |
            (np.shape(lac_tof)[1] != np.shape(calibration_matrix)[1])):
            print(
                '!!!!!!! WARNING CALIBRATION MATRIX HAS NOT SAME SIZE OF IMAGE !!!!!!!!!!!!!!'
            )

    if (any(debug_idx)):  #testing on a single pixel
        lac = lac_tof[debug_idx[0], debug_idx[1], :]

        if (calibration_matrix.any()):
            lambd = rt.tof2l_t0k(
                spectrum, calibration_matrix[debug_idx[0], debug_idx[1], 1],
                calibration_matrix[debug_idx[0], debug_idx[1], 0])
        else:
            lambd = spectrum

        WavelengthSelectiveRatio(lac,
                                 lambd,
                                 l1,
                                 l2,
                                 l1_w,
                                 l2_w,
                                 bool_SG,
                                 SG_w,
                                 SG_n,
                                 bool_print=True)
        return

    phase_ratio = np.zeros(np.shape(mymask))
    for i in tqdm(range(0, np.shape(mymask)[0])):
        for j in range(0, np.shape(mymask)[1]):
            if (mymask[i, j]):
                lac = lac_tof[i, j, :]

                if (calibration_matrix.any()):
                    lambd = rt.tof2l_t0k(spectrum, calibration_matrix[i, j, 1],
                                         calibration_matrix[i, j, 0])
                else:
                    lambd = spectrum

                try:
                    phi_fit = WavelengthSelectiveRatio(lac,
                                                       lambd,
                                                       l1,
                                                       l2,
                                                       l1_w,
                                                       l2_w,
                                                       bool_SG,
                                                       SG_w,
                                                       SG_n,
                                                       bool_print=False)
                    phase_ratio[i, j] = phi_fit['WSR']
                except:
                    print("Unexpected error at :", i, j)
                    phase_ratio[i, j] = -2.0

    if (bool_print):
        plt.figure()
        plt.imshow(phase_ratio, cmap='jet'), plt.title('Phase ratio')
        plt.colorbar()
        plt.show(), plt.close()
    if (bool_save):
        np.save('phase_ratio.npy', phase_ratio)
    return {'phase_ratio': phase_ratio}