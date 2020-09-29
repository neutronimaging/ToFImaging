[Return to table of contents](index.md)<br/>
This package contains python functions for the fitting of a single Bragg edge pattern.
Currently implemented are:
    - Advanced Bragg Edge Fitting: 7 parameters fit as in Ramadhan et al. [2019] {a_0,b_0,a_{hkl},b_{hkl},λ_{hkl},σ,τ}
    ![img](https://latex.codecogs.com/svg.latex?T(\lambda)=\exp{[-a_0+b_0\lambda]}\times(\exp{[-(a_{hkl}+b_{hkl}\lambda)]+\{1-\exp{[-(a_{hkl}+b_{hkl}\lambda)]\})\times\frac{1}{2}[\mathrm{erfc}(-\frac{\lambda-\lambda_{hkl}}{2^{1/2}\sigma})-\exp{(-\frac{\lambda-\lambda_{hkl}}{\tau}+\frac{\sigma^2}{2\tau^2})}\times\mathrm{erfc}(-\frac{\lambda-\lambda_{hkl}}{2^{1/2}\sigma}+\frac{\sigma}{\tau})])
    - Gaussian Bragg Edge Fitting: This methods takes the attenuation or transmission derivative and fits the edge to a Gaussian, returning centroid, height, and width of the Gaussian fit. The height is converted to the bragg edge slope and the bragg edge height is calculated as the integral of the gaussian fit.

    - edgefitting_1D: functions for edge fitting of 1D-arrays
        - AdvancedBraggEdgeFitting: fit two phases
            - INPUTS:
                - signal = ndarray of the spectrum containing the Bragg edge(s)
                - spectrum = spectrum, length of this ndarray must correspond to size of Tspectrum (lambda)
                - spectrum_range = range corresponding to lambda where to perform the fitting
                - est_pos = estimated bragg edge position (in spectrum_range dimension)
                - est_sigma = expected Gaussian broadening
                - est_alpha = expected decay constant (moderator property)
                - bool_smooth = set to True to perform Savitzky-Golay filtering of the transmission derivative
                - smooth_w = window size of S-G filter
                - smooth_n = order of S-G filter
                - bool_linear = flag to activate linear spectrum assumptions at the sides of the edge (otherwise exponential)
                - bool_print = flag to activate printing of figures

            - OUTPUTS: dictionary with the following fit in the dimension of the mask
                - 't0' = fitted bragg edge position
                - 'sigma' = fitted Gaussian broadening
                - 'alpha' = fitted decay constant (moderator property)
                - 'a1, a2, a5, a6' = parameters for spectrum besides the edge: a1 and a2 before, a5 and a6 after the edge
                - 'final_result' = fitting result after 7th iteration
                - 'fitted_data' = final fitted spectrum 
                - 'pos_extrema' = extrema of the bragg edges
                - 'height' = fitted height of the bragg edge

        - GaussianBraggEdgeFitting: fit three phases
            - INPUTS:
                - lac = 1d array the attenuation -log(I/I0) (1darray) [REQUIRED]
                - spectrum = spectrum range corresponding to lac (1darray) [REQUIRED]
                - phase1lac = lac of the phase 1 (1darray) [REQUIRED]
                - phase2lac = lac of the phase 2 (1darray) [REQUIRED]
                - phase3lac = lac of the phase 3 (1darray) [REQUIRED]
                - spectrum_phase = spectrum range corresponding to phase1lac,phase2lac and phase3lac (1darray) [REQUIRED]
                - lambda_range_norm = lambda range where to normalize spectra ([lambda1, lambda2]) [REQUIRED]
                - lambda_range_edges = lambda range where to do the fitting ([lambda1, lambda2]) [REQUIRED]
                - est_f1 = estimate phase 1 weight [Default = 0.333]
                - est_f2 = estimate phase 2 weight [Default = 0.333]
                - est_f3 = estimate phase 3 weight [Default = 0.334]
                - method = fitting method [Default = 'least_squares']
                - bool_SG = set to True to perform Savitzky-Golay filtering [Default = False]
                - SG_w = window size of S-G filter [Default = 5]
                - SG_n = order of S-G filter [Default = 1]
                - bool_print = set to True to print output [Default = False]

            - OUTPUTS: dictionary with the following fit in the dimension of the mask
                - 'phi1' : phase 1 weight
                - 'phi2' : phase 2 weight
                - 'phi3' : phase 3 weight

    - edgefitting_2D: functions for phase fitting of 2D stack of TOF data in the form of 3darray (x,y,lambda)
        - GaussianBraggEdgeFitting_2D: fit two phases
            - INPUTS:
                - lac_tof = 3darray the attenuation -log(I/I0) (x,y,lambda) [REQUIRED]
                - spectrum = spectrum range corresponding to lac (1darray) [REQUIRED]
                - phase1lac = lac of the phase 1 (1darray) [REQUIRED]
                - phase2lac = lac of the phase 2 (1darray) [REQUIRED]
                - spectrum_phase = spectrum range corresponding to phase1lac and phase2lac (1darray) [REQUIRED]
                - lambda_range_norm = lambda range where to normalize spectra ([lambda1, lambda2]) [REQUIRED]
                - lambda_range_edges = lambda range where to do the fitting ([lambda1, lambda2]) [REQUIRED]
                - calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]). Will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y) [Default = np.ndarray([0])]
                - mask = mask of where to perform the fit (x,y) [Default = np.ndarray([0])]
                - auto_mask = if True and mask is not given will automatically mask the region based on the mask_thresh thresholds (bool) [Default = True]
                - mask_thresh = low and high threshold for the automatic mask ([thresh_low, thresh_high]) [Default = np.ndarray([0])]
                - est_phi = estimate phase 1 weight [Default = 0.5]
                - method = fitting method [Default = 'least_squares']
                - bool_SG = set to True to perform Savitzky-Golay filtering (bool) [Default = False]
                - SG_w = window size of S-G filter [Default = 5]
                - SG_n = order of S-G filter [Default = 1]
                - bool_save = set to True to save output (bool) [Default = False]
                - bool_print = set to True to print output [Default = False]
                - debug_idx = pixel coordinates where to test the single pixel fitting ([pixel_x, pixel_y]) [Default = []]

            - OUTPUTS: dictionary with the following fit in the dimension of the mask
                - 'phase_ratio' : phase 1 weight

        - AdvancedBraggEdgeFitting_2D: fit three phases
            - INPUTS:
                - lac_tof = 3darray the attenuation -log(I/I0) (x,y,lambda) [REQUIRED]
                - spectrum = spectrum range corresponding to lac (1darray) [REQUIRED]
                - phase1lac = lac of the phase 1 (1darray) [REQUIRED]
                - phase2lac = lac of the phase 2 (1darray) [REQUIRED]
                - phase3lac = lac of the phase 3 (1darray) [REQUIRED]
                - spectrum_phase = spectrum range corresponding to phase1lac,phase2lac and phase3lac (1darray) [REQUIRED]
                - lambda_range_norm = lambda range where to normalize spectra ([lambda1, lambda2]) [REQUIRED]
                - lambda_range_edges = lambda range where to do the fitting ([lambda1, lambda2]) [REQUIRED]
                - calibration_matrix = calibration matrix with the coefficients to convert from spectrum to lambda size (x,y,[X0,k]). Will convert to lambda using formula Y = X0 + kX where X is spectrum for each pixel (x,y) [Default = np.ndarray([0])]
                - mask = mask of where to perform the fit (x,y) [Default = np.ndarray([0])]
                - auto_mask = if True and mask is not given will automatically mask the region based on the mask_thresh thresholds (bool) [Default = True]
                - mask_thresh = low and high threshold for the automatic mask ([thresh_low, thresh_high]) [Default = np.ndarray([0])]
                - est_f1 = estimate phase 1 weight [Default = 0.333]
                - est_f2 = estimate phase 2 weight [Default = 0.333]
                - est_f3 = estimate phase 3 weight [Default = 0.334]
                - method = fitting method [Default = 'least_squares']
                - bool_SG = set to True to perform Savitzky-Golay filtering (bool) [Default = False]
                - SG_w = window size of S-G filter [Default = 5]
                - SG_n = order of S-G filter [Default = 1]
                - bool_save = set to True to save output (bool) [Default = False]
                - bool_print = set to True to print output [Default = False]
                - debug_idx = pixel coordinates where to test the single pixel fitting ([pixel_x, pixel_y]) [Default = []]

            - OUTPUTS: dictionary with the following fit in the dimension of the mask
                - 'phase1_ratio' : phase 1 weight
                - 'phase2_ratio' : phase 2 weight
                - 'phase3_ratio' : phase 3 weight