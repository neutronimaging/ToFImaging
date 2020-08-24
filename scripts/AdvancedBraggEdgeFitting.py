import numpy as np
from numpy import pi, r_, math, random
import matplotlib.pyplot as plt
#from scipy import optimize
#from scipy.optimize import curve_fit
from scipy.special import erfc, erf
from lmfit import Model
#from lmfit.printfuncs import *
#import lmfit
#from numpy import loadtxt

#from scipy.signal import argrelextrema
from TOF_routines import find_nearest
#from TOF_routines import find_first
#from TOF_routines import find_last
#from TOF_routines import rotatedata
from TOF_routines import savitzky_golay

# def term0(t,a2,a6):
#     return  a2 * (t - a6)

# def term1(t,a2,a5,a6):
#     return ((a5 - a2) / 2) * (t - a6)

def term3(t,t0,sigma):
    return erfc(-((t-t0)/(sigma * math.sqrt(2))))

def term3_1(t,t0,sigma):
    return erf(-((t-t0)/(sigma * math.sqrt(2))))

def term4(t,t0,alpha,sigma):
    return np.exp(-((t-t0)/alpha) + ((sigma*sigma)/(2*alpha*alpha)))

def term5(t,t0,alpha,sigma):
    return erfc(-((t-t0)/(sigma * math.sqrt(2))) + sigma/alpha)

def term5_1(t,t0,alpha,sigma):
    return erf(-((t-t0)/(sigma * math.sqrt(2))) + sigma/alpha)

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

def BraggEdgeLinear_Attenuation(t,t0,alpha,sigma,a1,a2,a5,a6):
    return line_before(t,a5,a6)*B(t,t0,alpha,sigma)+line_after(t,a1,a2)*(1-B(t,t0,alpha,sigma))

def BraggEdgeExponential(t,t0,alpha,sigma,a1,a2,a5,a6):
    return exp_after(t,a1,a2) * ( exp_before(t,a5,a6)+ (1-exp_before(t,a5,a6)) * B(t,t0,alpha,sigma) )

def BraggEdgeExponential_Attenuation(t,t0,alpha,sigma,a1,a2,a5,a6):
    return exp_before(t,a5,a6) * ( exp_after(t,a1,a2)+ (1-exp_after(t,a1,a2)) * B(t,t0,alpha,sigma) )

def SG_filter(x, w=3, n=1):
    return savitzky_golay(x, w, n)

def GaussianBraggEdgeFitting(myspectrum, myTOF, myrange=0, est_pos=0, est_wid=0, est_h=0, bool_smooth=False, smooth_w = 3, smooth_n = 0, bool_print=False):
    
    def gaussian(x, amp, cen, wid):
        """1-d gaussian: gaussian(x, amp, cen, wid)"""
        return (amp / (np.sqrt(2*pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))
        
    if(myrange):
        myspectrum = myspectrum[myrange[0]:myrange[1]]
        myTOF = myTOF[myrange[0]:myrange[1]]
        
    if (bool_smooth):
        myspectrum = SG_filter(myspectrum,smooth_w,smooth_n)    
        
    d_signal = np.diff(myspectrum)
    dtof = myTOF[0:-1]
        
    if (bool_smooth):
        d_signal = SG_filter(d_signal,smooth_w,smooth_n)    
    
    if (bool_print):
        plt.figure()
        plt.plot(myTOF, myspectrum)
        plt.plot(dtof, d_signal)
    
    ## 1st Appoach
    #init_vals = [myspectrum[-2]-myspectrum[2], est_pos, 0.01]  # for [amp, cen, wid]
    #best_vals, covar = curve_fit(gaussian, dtof, d_signal)
    #fitted_data = gaussian(dtof,best_vals[0],best_vals[1],best_vals[2])        
    #t0 = best_vals[1]
    #edge_width = best_vals[2]

    ## 2nd Appoach
    method ='least_squares' # default and it implements the Levenberg-Marquardt
    gmodel = Model(gaussian)
    if(est_pos==0):
        est_pos = dtof[np.int(len(dtof)/2)]
    if(est_wid==0):
        est_wid = 0.01
    if(est_h==0):
        est_h = myspectrum[-2]-myspectrum[2]

    params = gmodel.make_params(cen=est_pos, amp=est_h, wid=est_wid)
    result = gmodel.fit(d_signal, params, x=dtof, method=method, nan_policy='propagate')
    t0 = result.best_values.get('cen')
    edge_width = result.best_values.get('wid')
    fitted_data = result.best_fit       
    
    id_low = find_nearest(dtof, t0-edge_width) # 3.7
    id_high = find_nearest(dtof, t0+edge_width) # 3.7
    edge_height = np.sum(fitted_data[id_low:id_high])
    edge_slope = result.best_values.get('amp')
    
    if (bool_print):
        print('t0 = ',t0)
        print('width = ',edge_width)
        print('id_low = ',id_low,'id_high = ',id_high,'height = ',edge_height)
        plt.plot(t0, myspectrum[find_nearest(dtof, t0)],'x', markeredgewidth=3, c='orange')
        plt.plot(t0-edge_width, myspectrum[find_nearest(dtof, t0-edge_width)],'+', markeredgewidth=3, c='orange')
        plt.plot(t0+edge_width, myspectrum[find_nearest(dtof, t0+edge_width)],'+', markeredgewidth=3, c='orange')
        plt.plot(dtof, fitted_data)
        plt.plot(t0, d_signal[find_nearest(dtof, t0)],'x', markeredgewidth=3, c='orange')
        plt.plot(t0-edge_width, d_signal[find_nearest(dtof, t0-edge_width)],'+', markeredgewidth=3, c='orange')
        plt.plot(t0+edge_width, d_signal[find_nearest(dtof, t0+edge_width)],'+', markeredgewidth=3, c='orange')
        plt.title('Bragg edge')
        plt.xlabel('Wavelenght [Å]')
        plt.ylabel('Transmission I/I$_{0}$')
        #plt.savefig('step1_fitting.pdf')
        plt.show()
        plt.close()
    return {'fitted_data':fitted_data, 't0':t0, 'edge_width':edge_width, 'edge_height':edge_height, 'edge_slope':edge_slope}

def AdvancedBraggEdgeFitting(myspectrum, myrange, myTOF, est_pos=0, est_sigma=1, est_alpha=1, bool_print=False, bool_average=False, bool_linear=False): 
## my range should be now the index position of the spectra that I want to study, est_pos is also the index position where the expected peak is

    ##INPUTS:
    #myspectrum: spectrum to fit
    #myrange: array of length 2 with index of the spectrum to fit    
    #myTOF: tof domain of the spectrum to fit
    #est_pos: expected bragg edge position (INDEX)
    #est_sigma: expected Gaussian broadening
    #est_alpha: expected decay constant (moderator property)
    #bool_print: flag to activate printing of figures
    #bool_average: flag to activate moving average across the spectrum (actualy SG filter)
    #bool_linear: flag to activate linear spectrum assumptions at the sides of the edge (otherwise exponential)
    #bool_transmission: flag on if specturm is transmission, if it is off it assumes it is absorption and converts to Transmission [obsolete]
    
    ##OUTPUTS:
    #t0: fitted bragg edge position
    #sigma: fitted Gaussian broadening
    #alpha: fitted decay constant (moderator property)
    #a1, a2, a5, a6: parameters for spectrum besides the edge: a1 and a2 before, a5 and a6 after the edge
    #final_result: fitting result after 7th iteration
    #fitted_data: final fitted spectrum 
    #pos_extrema: extrema of the bragg edges
    #height: fitted height of the bragg edge
    
    #get the part of the spectrum that I want to fit
    mybragg= myspectrum[myrange[0]:myrange[1]]
    t = myTOF[myrange[0]:myrange[1]]
    
    if(bool_average):
        mybragg = SG_filter(mybragg,3,0)
    
    if(est_pos==0):
        est_pos = np.argmax(SG_filter(np.diff(mybragg)))
    #else:
        #est_pos = est_pos-myrange[0] # I move the estimated position relative to the studied range, this is an index
   
    t0_f=t[est_pos] # this is the actual estimated first position in TOF [s]

    if (bool_print):
        plt.figure()
        plt.plot(t, mybragg)
        plt.plot(t0_f, mybragg[est_pos],'x', markeredgewidth=3, c='orange')
        plt.title('Bragg edge')
        plt.xlabel('Wavelenght [Å]')
        plt.ylabel('Transmission I/I$_{0}$')
        #plt.savefig('step1_fitting.pdf')
    
    t_before= t[0:est_pos]
    bragg_before=mybragg[0:est_pos]
    #     t_after= t[est_pos+int(est_pos*0.2):-1]
    #     bragg_after=mybragg[est_pos+int(est_pos*0.2):-1]
    t_after= t[est_pos+int(est_pos*0.2):-1]
    bragg_after=mybragg[est_pos+int(est_pos*0.2):-1]
    
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
            plt.plot(t,mybragg)
            plt.plot(t,interception_before+slope_before*t,'g')
            plt.plot(t,interception_after+slope_after*t,'r')
            plt.title('linear fitting before and after the given edge position')
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
            plt.plot(t,mybragg)
            plt.plot(t,interception_before+slope_before*t,'--r', label='fitted line before')
            plt.plot(t,interception_after+slope_after*t,'--g', label='fitted line after')
            plt.plot(t,exp_after(t,a1_f,a2_f),'g', label='fitted exp before')
            plt.plot(t,exp_combined(t,a1_f,a2_f,a5_f,a6_f),'r', label='fitted exp after')
            plt.xlabel('Wavelenght [Å]')
            plt.ylabel('Transmission I/I$_{0}$')
            plt.title('fitting before and after the given edge position')
            plt.legend()
            #plt.savefig('step2_fitting_legend.pdf')
            plt.plot(t, BraggEdgeExponential(t,t0_f,est_alpha,est_sigma,a1_f,a2_f,a5_f,a6_f))

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
    first_guess = gmodel.eval(params, t=t)
    if (bool_print):
        plt.figure()
        plt.plot(t,mybragg,label='data')
        plt.plot(t,first_guess,'--b',label='initial model')
        plt.title('initial BE with given parameters')
        plt.xlabel('Wavelenght [Å]')
        plt.ylabel('Transmission I/I$_{0}$')
        plt.legend()
        #plt.savefig('step3_fitting_legend.pdf')
    
    params['alpha'].vary = False
    params['sigma'].vary = False
    params['t0'].vary = False
    params['a2'].vary= False
    params['a5'].vary = False
    
    result1 = gmodel.fit(mybragg, params, t=t, method=method, nan_policy='propagate') 
    #    print(result1.fit_report())
    a1_f=result1.best_values.get('a1')
    a6_f=result1.best_values.get('a6')
    
    # 2nd round to fit a2 and a5 (slope of before and intercept of after)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['alpha'].vary = False
    params['sigma'].vary = False
    params['t0'].vary = False
    params['a1'].vary= False
    params['a6'].vary = False
    
    result2 = gmodel.fit(mybragg, params, t=t, method=method, nan_policy='propagate')
    a2_f = result2.best_values.get('a2')
    a5_f = result2.best_values.get('a5')
    
    # 3rd round to fit t0 (position)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['a2'].vary = False
    params['a5'].vary = False
    params['a1'].vary= False
    params['a6'].vary = False
    params['sigma'].vary= False
    params['alpha'].vary = False
    params['t0'].min = myTOF[myrange[0]]
    params['t0'].max = myTOF[myrange[1]]
    
    result3=gmodel.fit(mybragg, params, t=t, method=method, nan_policy='propagate')
    t0_f=result3.best_values.get('t0')
    
    # 4th round to fit sigma, alpha and t0 (broadening, decay and position)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['a2'].vary = False
    params['a5'].vary = False
    params['a1'].vary= False
    params['a6'].vary = False
    params['t0'].min = myTOF[myrange[0]]
    params['t0'].max = myTOF[myrange[1]]
    params['alpha'].min = 0.0
    params['alpha'].max = 1.5
    params['sigma'].min = 0.0
    params['sigma'].max = 1.5
    
    result4=gmodel.fit(mybragg, params, t=t, nan_policy='propagate',method=method)
    sigma_f=result4.best_values.get('sigma')
    alpha_f=result4.best_values.get('alpha')
    t0_f=result4.best_values.get('t0')
    
    # 5th round to fit a1 a2 a5 and a6 (slope of intercepts of before and after)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['t0'].vary = False
    params['sigma'].vary = False
    params['alpha'].vary= False
   
    result5 = gmodel.fit(mybragg, params, t=t, nan_policy='propagate', method=method)
    a1_f =result5.best_values.get('a1')
    a2_f = result5.best_values.get('a2')
    a5_f = result5.best_values.get('a5')
    a6_f = result5.best_values.get('a6')

    # 6th round to fit sigma alpha and t0 (broadening, decay and position)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['a2'].vary = False
    params['a5'].vary = False
    params['a1'].vary= False
    params['a6'].vary = False
    params['t0'].min = myTOF[myrange[0]]
    params['t0'].max = myTOF[myrange[1]]
    params['alpha'].min = 0.0
    params['alpha'].max =1.5
    params['sigma'].min = 0.0
    params['sigma'].max =1.5
    
    result6= gmodel.fit(mybragg, params, t=t, nan_policy='propagate', method=method)   
    sigma_f=result6.best_values.get('sigma')
    alpha_f=result6.best_values.get('alpha')
    t0_f=result6.best_values.get('t0')
    
    # 7th round to fit sigma alpha and t0 a1 a2 a5 a6 (all)
    params = gmodel.make_params(t0=t0_f,sigma=sigma_f, alpha=alpha_f, a1=a1_f, a2=a2_f, a5=a5_f, a6=a6_f)
    params['t0'].min = myTOF[myrange[0]]
    params['t0'].max = myTOF[myrange[1]]
    params['alpha'].min = 0.0
    params['alpha'].max =1.5
    params['sigma'].min = 0.0
    params['sigma'].max =1.5

    result7 = gmodel.fit(mybragg, params, t=t, nan_policy='propagate', method=method)    
    t0_f=result7.best_values.get('t0')
    sigma_f=result7.best_values.get('sigma')
    alpha_f=result7.best_values.get('alpha')
    a1_f =result7.best_values.get('a1')
    a2_f = result7.best_values.get('a2')
    a5_f = result7.best_values.get('a5')
    a6_f = result7.best_values.get('a6')
    
    if (bool_print):
        print(result7.fit_report())
        print(result7.covar)    
        print('bool value, Boolean for whether error bars were estimated by fit.', result7.errorbars)
        print(result7.ci_out) # print out the interval confidence
   
    #    Get the extrema for edge height fitting
    fitted_data = result7.best_fit
    pos_extrema = []
    
    ## Attempt n.1 -------- Here I was searching the last 0 value and the first value with 1.0 in the step function, however is it not a robust solution
    #     step_function = B(t,t0_f,alpha_f,sigma_f)
    #     min_pos = find_last(step_function,0.0)
    #     pos_extrema.append(min_pos)
    #     max_pos = find_first(step_function,0.99)
    #     pos_extrema.append(max_pos)
    
    if (bool_linear):
        fit_before = line_before(t,a5_f,a6_f)
        fit_after = line_after(t,a1_f,a2_f)
    else:
        fit_before = exp_combined(t,a1_f,a2_f,a5_f,a6_f)
        fit_after = exp_after(t,a1_f,a2_f)

    #fit_edge = B(t,t0_f,alpha_f,sigma_f)
    
    if (bool_print):
        plt.figure()
        plt.plot(t,fit_before,'o-')
        plt.plot(t,fit_after,'o-')
        plt.plot(t,fitted_data,'.-')
    
    index_t0 = find_nearest(t,t0_f)

    ## Attempt n.2 ------This is Florencia's approach and most used in the community: be careful that sometimes this gives an overestimation depending on the slope of the lines before and after
    #     pos_extrema.append(fit_before[index_t0])
    #     pos_extrema.append(fit_after[index_t0])
    #     height = np.abs(fit_after[index_t0]-fit_before[index_t0])

    ## Attempt n.3 ------ This approach is based on the difference between the fit before and after the edge and the fitted data itself. so far, it gives the nicest results on the calibration sample, however the value used as threshold is not general and should probably be adjusted from case to case. So again it is not yet the final solution

    for i in range(0, len(fitted_data)):
#         print(i,(fitted_data[i]-fit_before[i]))
        if (np.abs(fitted_data[i]-fit_before[i])>1e-4):            
            pos_extrema.append(i-1)
            break

    for i in range(len(fitted_data)-1,0,-1): # here I am moving backwards
#         print(i,(fitted_data[i]-fit_after[i]))
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

    height = np.abs(mybragg[pos_extrema[0]]-mybragg[pos_extrema[1]])  
    
    if (bool_print):
#         print('first iteration: ' ,result1.fit_report())
#         print('second iteration: ', result2.fit_report())
#         print('third iteration: ', result3.fit_report())
#         print('fourth iteration: ', result4.fit_report())
#         print('fifth iteration: ', result5.fit_report())
#         print('sixth iteration: ', result6.fit_report())
        print(pos_extrema)
        plt.figure()
        plt.plot(t, mybragg)
        #     plt.plot(t, result7.init_fit, 'k--')

        plt.plot(t, result1.best_fit, '--', color='gray', label='intermediate steps')
        plt.plot(t, result1.init_fit, '--', color='gray')
        plt.plot(t, result2.best_fit, '--', color='gray')
        plt.plot(t, result3.best_fit, '--', color='gray')
        plt.plot(t, result4.best_fit, '--', color='gray')
        plt.plot(t, result5.best_fit, '--', color='gray')
        plt.plot(t, result6.best_fit, '--', color='gray')
        plt.plot(t, result7.best_fit, 'r', linewidth='1.5', label='final fit')
        plt.legend()
        plt.xlabel('Wavelenght [Å]')
        plt.ylabel('Transmission I/I$_{0}$')

        plt.plot(t0_f,result7.best_fit[index_t0],'ok')
        plt.plot(t[pos_extrema[0]],result7.best_fit[pos_extrema[0]],'ok')
        plt.plot(t[pos_extrema[1]],result7.best_fit[pos_extrema[1]],'ok')
        plt.title('edge fitting and estimated edge position')
        plt.show()
    
    return {'t0':t0_f, 'sigma':sigma_f, 'alpha':alpha_f, 'a1':a1_f, 'a2':a2_f,'a5':a5_f, 'a6':a6_f, 'final_result':result7, 'fitted_data':fitted_data, 'pos_extrema':pos_extrema, 'height':height}

def phase_fraction_fitting(lac,spectrum_lambda,phase1lac,phase2lac,phase_spectrum,lambda_range_norm,lambda_range_edges,bool_plot=0): 
    ##INPUTS:
    #lac: measured LAC spectrum to fit
    #spectrum_lambda: lambda spectrum of the lac
    #phase1lac: LAC of the 1st phase to fit
    #phase2lac: LAC of the 2nd phase to fit
    #phase_spectrum: lambda spectrum of the theoretical phases LAC
    #lambda_range_norm: lambda range used to normalize data
    #lambda_range_edges: lambda range where to perform the fitting
    #bool_plot: flag to activate plotting
    
    ##OUTPUTS:
    #phi: fitted phase fraction

    #cutting to selected range for normalization
    idx_low = find_nearest(spectrum_lambda,lambda_range_norm[0])
    idx_high = find_nearest(spectrum_lambda,lambda_range_norm[1])
    lac = lac[idx_low:idx_high]
    spectrum_lambda = spectrum_lambda[idx_low:idx_high]

    idx_low = find_nearest(phase_spectrum,lambda_range_norm(1))
    idx_high = find_nearest(phase_spectrum,lambda_range_norm(2))
    lambda_table = phase_spectrum[idx_low:idx_high]
    phase1lac = phase1lac[idx_low:idx_high]
    phase2lac = phase2lac[idx_low:idx_high]

    phase1lac = np.interp(lambda_table,phase1lac,spectrum_lambda)
    phase2lac =  np.interp(lambda_table,phase2lac,spectrum_lambda)
    phase1lac = phase1lac*(np.nanmean(lac)/np.nanmean(phase1lac))
    phase2lac = phase2lac*(np.nanmean(lac)/np.nanmean(phase2lac))

    # initial conditions
    idx_l = find_nearest(spectrum_lambda,lambda_range_edges[0])
    idx_h = find_nearest(spectrum_lambda,lambda_range_edges[1])

    # method='trust_exact'
    # method='nelder' #not bad
    # method='differential_evolution' # needs bounds
    # method='basinhopping' # not bad
    # method='lmsquare' # this should implement the Levemberq-Marquardt but it says Nelder-Mead method (which should be Amoeba)
    method ='least_squares' # default and it implements the Levenberg-Marquardt
    def phase_linearcomb(ph1,ph2,f):
        return f*ph1+(1-f)*ph2

    gmodel = Model(phase_linearcomb)
    params = gmodel.make_params(f = 0.5, ph1 = phase1lac[idx_l:idx_h], ph2 = phase2lac[idx_l:idx_h])    
    params['ph1'].vary = False
    params['ph2'].vary = False
    params['f'].vary = True
    params['f'].min = 0.0
    params['f'].max = 1.0
    result = gmodel.fit(data = lac[idx_l:idx_h], params = params, method=method, nan_policy='propagate')
    phi=result.best_values.get('f')    

    return {'phi': phi}