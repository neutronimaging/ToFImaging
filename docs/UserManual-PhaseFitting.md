[Return to table of contents](index.md)<br/>
This package contains python functions for the fitting of phase transformation of Bragg attenuation patterns.
This assumes the measured linear attenuation coefficient μ, is a linear combination of a phases with a given a-priori reference linear attenuation coefficients φ

![img](https://latex.codecogs.com/svg.latex?\mu%3D\sum_i%20f_i%20\phi_i)

Currently implemented are the solvers for the 2 phases case and 3 phases case.

    - phasefitting_1D: functions for phase fitting of 1D-arrays
        - phase_ratio_linearcomb:
            - INPUTS:
                - lac = 1d array the attenuation -log(I/I0) TOF images (x,y,lambda) 
                - spectrum = spectrum, length of this ndarray must correspond to size of lac_tof(lambda)
                - phase1lac = lac of the phase 1
                - phase2lac = lac of the phase 2
                - spectrum_phase = spectrum corresponding to phase 1 and 2
                - lambda_range_norm = lambda range where to normalize spectra
                - lambda_range_edges = lambda range where to do the fitting
                - est_phi = estimate phase 1 weight
                - method = fitting method
                - bool_SG = set to True to perform Savitzky-Golay filtering of the transmission derivative
                - SG_w = window size of S-G filter
                - SG_n = order of S-G filter
                - bool_print = set to True to print output

            - OUTPUTS: dictionary with the following fit in the dimension of the mask
                - 'phi' : phase1 weight
                
        - phase_ratio_linearcomb_three:

    - phasefitting_2D: functions for phase fitting of 1D-arrays
        - phase_ratio_linearcomb_2D:
        - phase_ratio_linearcomb_2D_Calib_matrix:
        - phase_ratio_linearcomb_three:
        - phase_ratio_linearcomb_three_2D_Calib_matrix: