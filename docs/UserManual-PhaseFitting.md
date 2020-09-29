[Return to table of contents](index.md)<br/>
This package contains python functions for the fitting of phase transformation of Bragg attenuation patterns.
This assumes the measured linear attenuation coefficient μ, is a linear combination of a phases with a given a-priori reference linear attenuation coefficients φ

![img](https://latex.codecogs.com/svg.latex?\mu%3D\sum_i%20f_i%20\phi_i)

Currently implemented are the solvers for the 2 phases case (i=1,2) and 3 phases case (i=1,2,3).

    - phasefitting_1D: functions for phase fitting of 1D-arrays
        - phase_ratio_linearcomb:
            - INPUTS:
                - lac = 1d array the attenuation -log(I/I0) TOF images (lambda) [REQUIRED]
                - spectrum = spectrum, length of this ndarray must correspond to size of lac_tof (lambda) [REQUIRED]
                - phase1lac = lac of the phase 1 (lambda) [REQUIRED]
                - phase2lac = lac of the phase 2 (lambda) [REQUIRED]
                - spectrum_phase = spectrum corresponding to phase 1 and 2 (lambda) [REQUIRED]
                - lambda_range_norm = lambda range where to normalize spectra ([lambda1, lambda2]) [REQUIRED]
                - lambda_range_edges = lambda range where to do the fitting ([lambda1, lambda2]) [REQUIRED]
                - est_phi = estimate phase 1 weight [Default = 0.5]
                - method = fitting method [Default = 'least_squares']
                - bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative [Default = False]
                - SG_w = window size of S-G filter [Default = 5]
                - SG_n = order of S-G filter [Default = 1]
                - bool_print = set to True to print output [Default = False]

            - OUTPUTS: dictionary with the following fit in the dimension of the mask
                - 'phi' : phase1 weight

        - phase_ratio_linearcomb_three:
            - INPUTS:
                - lac = 1d array the attenuation -log(I/I0) TOF images (lambda) [REQUIRED]
                - spectrum = spectrum, length of this ndarray must correspond to size of lac_tof (lambda) [REQUIRED]
                - phase1lac = lac of the phase 1 (lambda) [REQUIRED]
                - phase2lac = lac of the phase 2 (lambda) [REQUIRED]
                - phase3lac = lac of the phase 3 (lambda) [REQUIRED]
                - spectrum_phase = spectrum corresponding to phase 1 and 2 ([lambda1, lambda2]) [REQUIRED]
                - lambda_range_norm = lambda range where to normalize spectra ([lambda1, lambda2]) [REQUIRED]
                - lambda_range_edges = lambda range where to do the fitting ([lambda1, lambda2]) [REQUIRED]
                - est_f1 = estimate phase 1 weight [Default = 0.333]
                - est_f2 = estimate phase 2 weight [Default = 0.333]
                - est_f3 = estimate phase 3 weight [Default = 0.334]
                - method = fitting method [Default = 'least_squares']
                - bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative [Default = False]
                - SG_w = window size of S-G filter [Default = 5]
                - SG_n = order of S-G filter [Default = 1]
                - bool_print = set to True to print output [Default = False]

            - OUTPUTS: dictionary with the following fit in the dimension of the mask
                - 'phi1' : phase 1 weight
                - 'phi2' : phase 2 weight
                - 'phi3' : phase 3 weight

    - phasefitting_2D: functions for phase fitting of 1D-arrays
        - phase_ratio_linearcomb_2D:
        - phase_ratio_linearcomb_2D_Calib_matrix:
        - phase_ratio_linearcomb_three:
        - phase_ratio_linearcomb_three_2D_Calib_matrix: