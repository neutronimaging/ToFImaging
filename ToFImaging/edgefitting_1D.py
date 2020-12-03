import numpy as np
from numpy import pi, r_, math, random
import matplotlib.pyplot as plt
from scipy import special
from lmfit import Model

from ToFImaging.reduction_tools import find_nearest
from ToFImaging.reduction_tools import savitzky_golay as SG_filter

def AdvancedBraggEdgeFitting(signal,spectrum,spectrum_range=[],est_pos=0,est_sigma=1,est_alpha=1,bool_smooth=False,smooth_w=5,smooth_n=1,bool_linear=False,bool_print=False): 
    """ Performs Bragg edge fitting with gaussian model to an ndarray containing the signal with the length of the spectrum (could be lambda, tof or bin index)

    INPUTS:
    signal = ndarray of the spectrum containing the Bragg edge(s)
    spectrum = spectrum, length of this ndarray must correspond to size of Tspectrum (lambda)
    spectrum_range = range corresponding to lambda where to perform the fitting
    est_pos = estimated bragg edge position (in spectrum_range dimension)
    est_sigma = expected Gaussian broadening
    est_alpha = expected decay constant (moderator property)
    bool_smooth = set to True to perform Savitzky-Golay filtering of the transmission derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    bool_linear = flag to activate linear spectrum assumptions at the sides of the edge (otherwise exponential)
    bool_print = flag to activate printing of figures
    
    OUTPUTS:
    dictionary with the following fits
    't0' = fitted bragg edge position
    'sigma' = fitted Gaussian broadening
    'alpha' = fitted decay constant (moderator property)
    'a1, a2, a5, a6' = parameters for spectrum besides the edge: a1 and a2 before, a5 and a6 after the edge
    'final_result' = fitting result after 7th iteration
    'fitted_data' = final fitted spectrum 
    'pos_extrema' = extrema of the bragg edges
    'height' = fitted height of the bragg edge
    """     

    #----------------- FITTING FUNCTIONS---------------#
    # def term0(t,a2,a6):
    #     return  a2 * (t - a6)

    # def term1(t,a2,a5,a6):
    #     return ((a5 - a2) / 2) * (t - a6)

    # def BraggEdgeExponential_Attenuation(t,t0,alpha,sigma,a1,a2,a5,a6):
    #     return exp_before(t,a5,a6) * ( exp_after(t,a1,a2)+ (1-exp_after(t,a1,a2)) * B(t,t0,alpha,sigma) )

    # def term3_1(t,t0,sigma):
    #     return math.erf(-((t-t0)/(sigma * math.sqrt(2))))

    # def term5_1(t,t0,alpha,sigma):
    #     return math.erf(-((t-t0)/(sigma * math.sqrt(2))) + sigma/alpha)

    # def BraggEdgeLinear_Attenuation(t,t0,alpha,sigma,a1,a2,a5,a6):
    #     return line_before(t,a5,a6)*B(t,t0,alpha,sigma)+line_after(t,a1,a2)*(1-B(t,t0,alpha,sigma))

    def term3(t,t0,sigma):
        return special.erfc(-((t-t0)/(sigma * math.sqrt(2))))
        
    def term4(t,t0,alpha,sigma):
        return np.exp(-((t-t0)/alpha) + ((sigma*sigma)/(2*alpha*alpha)))

    def term5(t,t0,alpha,sigma):
        return special.erfc(-((t-t0)/(sigma * math.sqrt(2))) + sigma/alpha)

    def line_after(t,a1,a2):
        return a1+a2*t

    def line_before(t,a5,a6):
        return a5+a6*t

    def exp_after(t,a1,a2):
        return np.exp(-(a1+a2*t))

    def exp_before(t,a5,a6):
        return np.exp(-(a5+a6*t))

    def exp_combined(t,a1,a2,a5,a6):
        return exp_after(t,a1,a2)*exp_before(t,a5,a6)

    def B(t,t0,alpha,sigma):
        edge = 0.5*(term3(t,t0,sigma) - term4(t,t0,alpha,sigma)* term5(t,t0,alpha,sigma))
        return (edge)

    def BraggEdgeLinear(t,t0,alpha,sigma,a1,a2,a5,a6):
        return line_after(t,a1,a2)*B(t,t0,alpha,sigma)+line_before(t,a5,a6)*(1-B(t,t0,alpha,sigma))

    def BraggEdgeExponential(t,t0,alpha,sigma,a1,a2,a5,a6):
        return exp_after(t,a1,a2) * ( exp_before(t,a5,a6)+ (1-exp_before(t,a5,a6)) * B(t,t0,alpha,sigma) )

    #-----------------FITTING PART---------------#
    t=spectrum #renamed to t for convenience but could be lambda or bin index
    #get the part of the spectrum that I want to fit
    if(spectrum_range):
        idx_low = find_nearest(spectrum,spectrum_range[0])
        idx_high = find_nearest(spectrum,spectrum_range[1])        
        signal = signal[idx_low:idx_high]
        t = t[idx_low:idx_high]
    
    if(bool_smooth):
        signal = SG_filter(signal,smooth_w,smooth_n)
    
    if(est_pos==0):
        est_pos = np.argmax(SG_filter(np.diff(signal)))
    else:
        est_pos = find_nearest(t,est_pos)        
    t0_f=t[est_pos] # this is the actual estimated first position in TOF [s]

    if (bool_print):
        plt.figure()
        plt.plot(t, signal)
        plt.plot(t0_f, signal[est_pos],'x', markeredgewidth=3, c='orange')
        plt.title('Bragg edge pattern and initial guess'), plt.xlabel('Wavelenght [Å]'), plt.ylabel('Transmission I/I$_{0}$')
        plt.show()
        plt.close()
        #plt.savefig('step1_fitting.pdf')
    
    t_before= t[0:est_pos]
    bragg_before=signal[0:est_pos]
    t_after= t[est_pos+int(est_pos*0.2):-1]
    bragg_after=signal[est_pos+int(est_pos*0.2):-1]
    
    #first step: estimate the linear or exponential function before and after the Bragg Edge
    [slope_before, interception_before] = np.polyfit(t_before, bragg_before, 1)
    [slope_after, interception_after] = np.polyfit(t_after, bragg_after, 1)
    #first guess of paramters
    a2_f=slope_after
    a5_f=interception_before
    a6_f=slope_before
    a1_f=interception_after    
    if (bool_linear):
        gmodel = Model(BraggEdgeLinear)
        if (bool_print):
            plt.figure()
            plt.plot(t_before,bragg_before,'.g')
            plt.plot(t_after,bragg_after,'.r')
            plt.plot(t,signal)
            plt.plot(t,interception_before+slope_before*t,'g')
            plt.plot(t,interception_after+slope_after*t,'r')
            plt.title('Linear fitting before and after the given edge position'), plt.xlabel('Wavelenght [Å]'), plt.ylabel('Transmission I/I$_{0}$')
            plt.show()
            plt.close()
    else:
        exp_model_after = Model(exp_after)
        params = exp_model_after.make_params(a1=a1_f, a2=a2_f)
        result_exp_model_after = exp_model_after.fit(bragg_after,params,t=t_after)
        a1_f=result_exp_model_after.best_values.get('a1')
        a2_f=result_exp_model_after.best_values.get('a2')
        
        exp_model_before = Model(exp_combined)
        params = exp_model_before.make_params(a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
        params['a1'].vary = False
        params['a2'].vary = False
        result_exp_model_before = exp_model_before.fit(bragg_before,params,t=t_before)
        a5_f=result_exp_model_before.best_values.get('a5')
        a6_f=result_exp_model_before.best_values.get('a6')
        gmodel = Model(BraggEdgeExponential)
        
        if (bool_print):
            plt.figure()
            plt.plot(t_before,bragg_before,'.r', label ='int point')
            plt.plot(t_after,bragg_after,'.g', label='int point')
            plt.plot(t,signal)
            plt.plot(t,interception_before+slope_before*t,'--r', label='fitted line before')
            plt.plot(t,interception_after+slope_after*t,'--g', label='fitted line after')
            plt.plot(t,exp_after(t,a1_f,a2_f),'g', label='fitted exp before')
            plt.plot(t,exp_combined(t,a1_f,a2_f,a5_f,a6_f),'r', label='fitted exp after')
            plt.xlabel('Wavelenght [Å]'), plt.ylabel('Transmission I/I$_{0}$') 
            plt.title('Exponential and line fitting adjacent to Bragg edge')
            plt.legend()
            plt.plot(t, BraggEdgeExponential(t,t0_f,est_alpha,est_sigma,a1_f,a2_f,a5_f,a6_f))
            plt.show()
            plt.close()

    sigma_f = est_sigma
    alpha_f = est_alpha

    # method='trust_exact'
    # method='nelder' #not bad
    # method='differential_evolution' # needs bounds
    # method='basinhopping' # not bad
    # method='lmsquare' # this should implement the Levemberq-Marquardt but it says Nelder-Mead method (which should be Amoeba)
    method ='least_squares' # default and it implements the Levenberg-Marquardt
    
    # 1st round to fit a1 and a6 (intercept of before and slope of after)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    # first_guess = gmodel.eval(params, t=t)
    
    params['alpha'].vary = False
    params['sigma'].vary = False
    params['t0'].vary = False
    params['a2'].vary = False
    params['a5'].vary = False
    
    result1 = gmodel.fit(signal, params, t=t, method=method, nan_policy='propagate') 
    a1_f = result1.best_values.get('a1')
    a6_f = result1.best_values.get('a6')
    
    #2nd round to fit a2 and a5 (slope of before and intercept of after)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['alpha'].vary = False
    params['sigma'].vary = False
    params['t0'].vary = False
    params['a1'].vary = False
    params['a6'].vary = False
    
    result2 = gmodel.fit(signal, params, t=t, method=method, nan_policy='propagate')
    a2_f = result2.best_values.get('a2')
    a5_f = result2.best_values.get('a5')
    
    #3rd round to fit t0 (position)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['a2'].vary = False
    params['a5'].vary = False
    params['a1'].vary = False
    params['a6'].vary = False
    params['sigma'].vary = False
    params['alpha'].vary = False
    params['t0'].min = t[0]
    params['t0'].max = t[-1]
    
    result3 = gmodel.fit(signal, params, t=t, method=method, nan_policy='propagate')
    t0_f = result3.best_values.get('t0')
    
    # 4th round to fit sigma, alpha and t0 (broadening, decay and position)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['a2'].vary = False
    params['a5'].vary = False
    params['a1'].vary = False
    params['a6'].vary = False
    params['t0'].min = t[0]
    params['t0'].max = t[-1]
    params['alpha'].min = 0.0
    params['alpha'].max = 1.5
    params['sigma'].min = 0.0
    params['sigma'].max = 1.5
    
    result4 = gmodel.fit(signal, params, t=t, nan_policy='propagate',method=method)
    sigma_f = result4.best_values.get('sigma')
    alpha_f = result4.best_values.get('alpha')
    t0_f = result4.best_values.get('t0')
    
    #5th round to fit a1 a2 a5 and a6 (slope of intercepts of before and after)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['t0'].vary = False
    params['sigma'].vary = False
    params['alpha'].vary = False
   
    result5 = gmodel.fit(signal, params, t=t, nan_policy='propagate', method=method)
    a1_f = result5.best_values.get('a1')
    a2_f = result5.best_values.get('a2')
    a5_f = result5.best_values.get('a5')
    a6_f = result5.best_values.get('a6')

    #6th round to fit sigma alpha and t0 (broadening, decay and position)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['a2'].vary = False
    params['a5'].vary = False
    params['a1'].vary = False
    params['a6'].vary = False
    params['t0'].min = t[0]
    params['t0'].max = t[-1]
    params['alpha'].min = 0.0
    params['alpha'].max = 1.5
    params['sigma'].min = 0.0
    params['sigma'].max = 1.5
    
    result6 = gmodel.fit(signal, params, t=t, nan_policy='propagate', method=method)   
    sigma_f = result6.best_values.get('sigma')
    alpha_f = result6.best_values.get('alpha')
    t0_f = result6.best_values.get('t0')
    
    #7th round to fit sigma alpha and t0 a1 a2 a5 a6 (all)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['t0'].min = t[0]
    params['t0'].max = t[-1]
    params['alpha'].min = 0.0
    params['alpha'].max = 1.5
    params['sigma'].min = 0.0
    params['sigma'].max = 1.5

    result7 = gmodel.fit(signal, params, t=t, nan_policy='propagate', method=method)    
    t0_f = result7.best_values.get('t0')
    sigma_f = result7.best_values.get('sigma')
    alpha_f = result7.best_values.get('alpha')
    a1_f = result7.best_values.get('a1')
    a2_f = result7.best_values.get('a2')
    a5_f = result7.best_values.get('a5')
    a6_f = result7.best_values.get('a6')
    
    if (bool_print):
        print(result7.fit_report())
        print(result7.covar)    
        print('bool value, Boolean for whether error bars were estimated by fit.', result7.errorbars)
        print(result7.ci_out) # print out the interval confidence
   
    #Get the extrema for edge height fitting
    fitted_data = result7.best_fit
    pos_extrema = []

    if (bool_linear):
        fit_before = line_before(t,a5_f,a6_f)
        fit_after = line_after(t,a1_f,a2_f)
    else:
        fit_before = exp_combined(t,a1_f,a2_f,a5_f,a6_f)
        fit_after = exp_after(t,a1_f,a2_f)

    index_t0 = find_nearest(t,t0_f)

    # This plots the fitted edge and the sides
    # if (bool_print): 
    #     plt.figure()
    #     plt.plot(t,fit_before,'o-')
    #     plt.plot(t,fit_after,'o-')
    #     plt.plot(t,fitted_data,'.-')
    #     plt.title('Final Fit'), plt.xlabel('Wavelenght [Å]'), plt.ylabel('Transmission I/I$_{0}$')

    ## Attempt n.1 -------- Here I was searching the last 0 value and the first value with 1.0 in the step function, however is it not a robust solution
    #     step_function = B(t,t0_f,alpha_f,sigma_f)
    #     min_pos = find_last(step_function,0.0)
    #     pos_extrema.append(min_pos)
    #     max_pos = find_first(step_function,0.99)
    #     pos_extrema.append(max_pos)

    #fit_edge = B(t,t0_f,alpha_f,sigma_f)

    ## Attempt n.2 ------This is Florencia's approach and most used in the community: be careful that sometimes this gives an overestimation depending on the slope of the lines before and after
    #     pos_extrema.append(fit_before[index_t0])
    #     pos_extrema.append(fit_after[index_t0])
    #     height = np.abs(fit_after[index_t0]-fit_before[index_t0])

    ## Attempt n.3 ------ This approach is based on the difference between the fit before and after the edge and the fitted data itself. so far, it gives the nicest results on the calibration sample, however the value used as threshold is not general and should probably be adjusted from case to case. So again it is not yet the final solution
    for i in range(0, len(fitted_data)):
        #print(i,(fitted_data[i]-fit_before[i]))
        if (np.abs(fitted_data[i]-fit_before[i])>1e-4):            
            pos_extrema.append(i-1)
            break

    for i in range(len(fitted_data)-1,0,-1): # here I am moving backwards
        #print(i,(fitted_data[i]-fit_after[i]))
        if (np.abs(fitted_data[i]-fit_after[i])>1e-3):            
            pos_extrema.append(i)
            break

    ## Attempt n.4 -- max and min before and after the estimated edge position, for the calibration sample works fine
    #     range_min = t[0:index_t0]
    #     range_max= t[index_t0:-1]
    #     min_fit = np.min(fitted_data[0:index_t0])
    #     max_fit = np.max(fitted_data[index_t0:-1])
    #     pos_min = find_nearest(fitted_data[0:index_t0], min_fit)
    #     pos_max = index_t0+find_nearest(fitted_data[index_t0:-1], max_fit)

    ## For other attempts have a look in the SENJU branch

    height = np.abs(signal[pos_extrema[0]]-signal[pos_extrema[1]])  
    
    if (bool_print):
        #print('first iteration: ' ,result1.fit_report())
        #print('second iteration: ', result2.fit_report())
        #print('third iteration: ', result3.fit_report())
        #print('fourth iteration: ', result4.fit_report())
        #print('fifth iteration: ', result5.fit_report())
        #print('sixth iteration: ', result6.fit_report())
        print(pos_extrema)
        plt.figure()
        plt.plot(t, signal)

        plt.plot(t, result1.best_fit, '--', color='gray', label='intermediate steps')
        plt.plot(t, result1.init_fit, '--', color='gray')
        plt.plot(t, result2.best_fit, '--', color='gray')
        plt.plot(t, result3.best_fit, '--', color='gray')
        plt.plot(t, result4.best_fit, '--', color='gray')
        plt.plot(t, result5.best_fit, '--', color='gray')
        plt.plot(t, result6.best_fit, '--', color='gray')
        plt.plot(t, result7.best_fit, 'r', linewidth='1.5', label='final fit')
        plt.legend(), plt.xlabel('Wavelenght [Å]'), plt.ylabel('Transmission I/I$_{0}$')

        plt.plot(t0_f,result7.best_fit[index_t0],'ok')
        plt.plot(t[pos_extrema[0]],result7.best_fit[pos_extrema[0]],'ok')
        plt.plot(t[pos_extrema[1]],result7.best_fit[pos_extrema[1]],'ok')
        plt.title('Edge fitting with intermediate steps and estimated edge position')
        plt.show()
        plt.close()
    
    return {'t0':t0_f, 'sigma':sigma_f, 'alpha':alpha_f, 'a1':a1_f, 'a2':a2_f,'a5':a5_f, 'a6':a6_f, 'final_result':result7, 'fitted_data':fitted_data, 'pos_extrema':pos_extrema, 'height':height}

def GaussianBraggEdgeFitting(signal,spectrum,spectrum_range=[],est_pos=0,est_wid=0,est_h=0,pos_BC=0,wid_BC=0,h_BC=0,bool_log=False,bool_smooth=False,smooth_w=5,smooth_n=1,interp_factor=0,bool_print=False):
    """ Performs Bragg edge fitting with gaussian model to an ndarray containing the signal with the length of the spectrum (could be lambda, tof or bin index)
    INPUTS:
    signal = ndarray of the spectrum containing the Bragg edge(s)
    spectrum = spectrum, length of this ndarray must correspond to size of Tspectrum(lambda)
    spectrum_range = range corresponding to lambda where to perform the fitting
    est_pos = estimated bragg edge position (in spectrum dimension)
    est_wid = estimated bragg edge width (in spectrum dimension)
    est_h = estimated bragg edge height (in spectrum dimension)
    pos_BC = boundary conditions for the bragg edge position fit
    wid_BC = boundary conditions for the bragg edge width fit
    h_BC = boundary conditions for the bragg edge height fit
    bool_log = set to True to perform log norm and convert to attenuation
    bool_smooth = set to True to perform Savitzky-Golay filtering of the signal derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    interp_factor = if set, this is the factor by which the number of bins are multiplied by interpolation
    bool_print = set to True to print output

    OUTPUTS:
    dictionary with the following fits
    'edge_position' : edge position 
    'edge_height': edge height 
    'edge_width': edge width  
    'edge_slope': edge slope 
    'median_image': median Transmission image in the selected lambda range
    """     
    def gaussian(x, amp, cen, wid):
        """1-d gaussian: gaussian(x, amp, cen, wid)"""
        return (amp / (np.sqrt(2*pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))
        
    if(spectrum_range):
        idx_low = find_nearest(spectrum,spectrum_range[0])
        idx_high = find_nearest(spectrum,spectrum_range[1])        
        signal = signal[idx_low:idx_high]
        spectrum = spectrum[idx_low:idx_high]
        
    if(bool_smooth):
        signal = SG_filter(signal,smooth_w,smooth_n)    
        
    d_signal = np.diff(signal)
    if(bool_log):
        d_signal = -d_signal

    d_spectrum = spectrum[0:-1]
        
    if (bool_smooth):
        d_signal = SG_filter(d_signal,smooth_w,smooth_n)    

    if(interp_factor):
        spectrum_n = np.linspace(d_spectrum[0],d_spectrum[-1],np.round(interp_factor*np.shape(d_spectrum)[0]))
        d_signal = np.interp(spectrum_n,d_spectrum,d_signal)
        d_spectrum = spectrum_n

    ## 2nd Appoach
    method ='least_squares' # default and it implements the Levenberg-Marquardt
    gmodel = Model(gaussian)
    if(est_pos==0):
        est_pos = d_spectrum[np.int(len(d_spectrum)/2)]
    if(est_wid==0):
        est_wid = 0.01
    if(est_h==0):
        est_h = signal[-2]-signal[2]

    params = gmodel.make_params(cen=est_pos, amp=est_h, wid=est_wid)
    if(pos_BC):
        params['cen'].min = pos_BC[0]
        params['cen'].max = pos_BC[1]
    else:
        params['cen'].min = spectrum[0]
        params['cen'].max = spectrum[-1]
    if(h_BC):
        params['amp'].min = h_BC[0]
        params['amp'].max = h_BC[1]
    else:
        params['amp'].min = 0
        params['amp'].max = 1e2
    if(wid_BC):
        params['wid'].min = wid_BC[0]
        params['wid'].max = wid_BC[1]
    else:
        params['wid'].min = 0
        params['wid'].max = 1e2

    result = gmodel.fit(d_signal, params, x=d_spectrum, method=method, nan_policy='propagate')
    t0 = result.best_values.get('cen')
    edge_width = result.best_values.get('wid')
    fitted_data = result.best_fit       
    
    id_low = find_nearest(d_spectrum, t0-edge_width) # 3.7
    id_high = find_nearest(d_spectrum, t0+edge_width) # 3.7
    edge_height = np.sum(fitted_data[id_low:id_high])
    edge_slope = result.best_values.get('amp')
    
    if (bool_print):
        print('Edge position = ', t0)
        print('Edge height = ', edge_height)
        print('Edge width = ', edge_width)
        print('idx_low = ', id_low,
              'idx_high = ', id_high)

        plt.figure()
        plt.subplot(2,1,1), 
        plt.plot(spectrum, signal)
        plt.plot(t0, signal[find_nearest(spectrum, t0)],'x', markeredgewidth=3, c='orange')
        plt.plot(t0-edge_width, signal[find_nearest(spectrum, t0-edge_width)],'+', markeredgewidth=3, c='orange')
        plt.plot(t0+edge_width, signal[find_nearest(spectrum, t0+edge_width)],'+', markeredgewidth=3, c='orange')
        plt.title('Bragg pattern'), plt.xlabel('Wavelenght [Å]')
        plt.subplot(2,1,2), 
        plt.plot(d_spectrum, d_signal), plt.plot(d_spectrum, fitted_data)
        plt.plot(t0, d_signal[find_nearest(d_spectrum, t0)],'x', markeredgewidth=3, c='orange')
        plt.plot(t0-edge_width, d_signal[find_nearest(d_spectrum, t0-edge_width)],'+', markeredgewidth=3, c='orange')
        plt.plot(t0+edge_width, d_signal[find_nearest(d_spectrum, t0+edge_width)],'+', markeredgewidth=3, c='orange')
        plt.title('Bragg pattern derivative'), plt.xlabel('Wavelenght [Å]')
        plt.tight_layout()
        plt.show()
        plt.close()

    return {'fitted_data': fitted_data,
            'edge_position': t0,
            't0': t0,
            'edge_width': edge_width,
            'edge_height': edge_height,
            'edge_slope': edge_slope}

#------------------------------ MICROSTRUCTURE FITTING ------------------------#
def MarchDollase(A,R,l,l_hkl,Nbeta=50,bool_plotPole=False):
    #A : angle between the hkl plane and preferred orientation
    theta = np.arcsin(l/l_hkl)
    beta = np.linspace(0,2*pi,Nbeta)
    theta = np.matlib.repmat(theta, Nbeta,1)
    beta = np.transpose(np.matlib.repmat(beta,np.shape(theta)[1],1))
    B = np.cos(A)*np.sin(theta)+np.sin(A)*np.multiply(np.cos(theta),np.sin(beta))
    P = np.sum(np.power(R**2*B**2+(1-B**2)/R,-3/2),axis=0)/(Nbeta)
    P[np.isnan(P)]=1
    if(bool_plotPole):
        plt.figure
        plt.imshow(np.power(R**2*B**2+(1-B**2)/R,-3/2))
        plt.show(), plt.close()
        plt.plot(l,P)
        plt.show(), plt.close()
    return P  

def SabinePrimaryExtinction(S,l,l_hkl):
    N = np.shape(l)[0]
    theta = np.arcsin(l/l_hkl)
    x = S**2*(l)**2
    E = np.zeros((N))
    for i in range(0,N):
        E_B = 1/np.sqrt(1+x[i])
        if(x[i]>1):
            E_L = np.sqrt(2/(pi*x[i]))*(1-1/(8*x[i])-3/(128*x[i]**2)-15/(1024*x[i]**3))
        else:
            E_L = 1-x[i]/2+x[i]**2/4-5*x[i]**3/48+7*x[i]**4/192
        E[i] = E_L*np.cos(theta[i])**2+E_B*np.sin(theta[i])**2
    return E

def TextureFitting(signal,spectrum,ref,ref_spectrum,spectrum_range=[],l_hkl1=1,l_hkl2=0,bool_MD=False,est_A1=0,est_R1=1,est_A2=0,est_R2=1,Nbeta=50,est_S=0,bool_smooth=False,smooth_w=5,smooth_n=1,bool_print=False):
    """ Performs texture fitting for up to two lattice planes 
    INPUTS:
    signal = ndarray of the spectrum containing the Bragg edge(s)
    spectrum = spectrum, length of this ndarray must correspond to size of Tspectrum(lambda)
    ref = reference spectrum for untextured cross sections
    ref_specturm = specturm of the reference
    spectrum_range = range corresponding to lambda where to perform the fitting
    l_hkl1 = wavelength of the first bragg edge
    l_hkl2 = wavelength of the second bragg edge (optional)
    bool_MD = switch for March-Dollase type texture fitting
    est_A1 = initial guess for preferred orientation for first bragg edge
    est_R1 = initial guess for anisotropy for first bragg edge
    est_A2 = initial guess for preferred orientation for second bragg edge
    est_R2 = initial guess for anisotropy for second bragg edge
    Nbeta = sampling number of the Pole figure (this size highly influence the fitting speed)
    est_S = initiaul guess for crystallite size (if set to 0 is not included in the fitting model)
    bool_smooth = set to True to perform Savitzky-Golay filtering of the signal derivative
    smooth_w = window size of S-G filter
    smooth_n = order of S-G filter
    bool_print = set to True to print output

    OUTPUTS:
    dictionary with the following fits
    'A1': preferred orientation for first bragg edge
    'R1': anisotropy for first bragg edge 
    'A2': preferred orientation for second bragg edge
    'R2': anisotropy for second bragg edge
    'S': crystallite size
    """   

    if(spectrum_range):
        idx_low = find_nearest(spectrum,spectrum_range[0])
        idx_high = find_nearest(spectrum,spectrum_range[1])
        signal = signal[idx_low:idx_high]
        spectrum = spectrum[idx_low:idx_high]

    if(bool_smooth):
        signal = SG_filter(signal,smooth_w,smooth_n)

    norm_sp=np.nanmean(signal)        
    ref_int = np.interp(spectrum,ref_spectrum,ref,left=np.nan,right=np.nan)
    
    def TexturedReference(ref_int,A1,R1,A2,R2,S):
        f = ref_int
        if(bool_MD): #if MarchDollase is set on apply it
            P1 = MarchDollase(A1,R1,spectrum,l_hkl1,Nbeta = Nbeta)
            P1[np.isnan(P1)]=1
            P=P1
            if(l_hkl2):
                P2 = MarchDollase(A2,R2,spectrum,l_hkl2,Nbeta = Nbeta)            
                P2[np.isnan(P2)]=1
                P=(P1+P2)/2
            f = np.multiply(f,P)                
        
        if(est_S): #if Crystallite size is set on apply it
            E1 = SabinePrimaryExtinction(S,spectrum,l_hkl1)
            E1[np.isnan(E1)]=1
            E=E1
            if(l_hkl2):
                E2 = SabinePrimaryExtinction(S,spectrum,l_hkl2)
                E2[np.isnan(E2)]=1
                E=(E1+E2)/2
            f = np.multiply(f,E)
        else: #When Crystallite size is not included data has to be rescaled
            f = f*norm_sp/np.nanmean(f)

        return f
        
    gmodel = Model(TexturedReference,independent_vars=['ref_int'])

    # if(bool_MD):
    #     params = gmodel.make_params(A1=est_A1,R1=est_R1)   
    #     if(l_hkl2):
    #         params = gmodel.make_params(A2=est_A2,R2=est_R2) 
    # if(est_S):
    #     params = gmodel.make_params(S=est_S)
    params = gmodel.make_params(A1=est_A1,R1=est_R1,A2=est_A2,R2=est_R2,S=est_S)

    if(bool_MD):
        params['A1'].min = 0.0
        params['A1'].max = np.pi/2
        params['R1'].min = 0.0
        params['R1'].max = 10.0
        if(l_hkl2):
            params['A2'].min = 0.0
            params['A2'].max = np.pi/2
            params['R2'].min = 0.0
            params['R2'].max = 10.0

    if(est_S):
        params['S'].min = 0.0
        params['S'].max = 1000.0

    method = 'least_squares'
    result = gmodel.fit(signal, params, ref_int = ref_int, method=method, nan_policy='propagate')

    #untextured default
    A1_fit=0
    R1_fit=1
    A2_fit=0
    R2_fit=1
    S_fit=0
    if(bool_MD):
        A1_fit=result.best_values.get('A1')    
        R1_fit=result.best_values.get('R1')         
        if(l_hkl2):
            A2_fit=result.best_values.get('A2')    
            R2_fit=result.best_values.get('R2')     

    if(est_S):
        S_fit=result.best_values.get('S')     

    if(bool_print):
        plt.plot(spectrum,signal,label='data')
        plt.plot(spectrum,ref_int,label='reference')
        plt.plot(spectrum,TexturedReference(ref_int,est_A1,est_R1,est_A2,est_R2,est_S),'--',label='Initial guess')
        plt.plot(spectrum,TexturedReference(ref_int,A1_fit,R1_fit,A2_fit,R2_fit,S_fit),'--',label='Fit')
        plt.legend(),plt.show(),plt.close()

    return {'A1': A1_fit, 'R1': R1_fit, 'A2': A2_fit, 'R2': R2_fit, 'S': S_fit}        