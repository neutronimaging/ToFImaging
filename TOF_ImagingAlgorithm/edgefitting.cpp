#include "edgefitting.h"
#include <lmcurve.h>
#include <QDebug>
#include <math/gradient.h>
#include <base/KiplException.h>

namespace ToFImagingAlgorithms {

/// class constructor, inputs:
/// n = number of parameters to be fitted, 7 for the complete lineshape, 3 for the simpliefies (gaussian of the signal gradient)
/// ef = lineshape type
edgefitting::edgefitting(int n, ToFImagingAlgorithms::eEdgeFunction ef)
{
    m_Npars = n;
    myfun = ef;
}

/// class deconstructor
edgefitting::~edgefitting()
{

}


/// initialize parameters by copying the inputs in the m_pars
void edgefitting::intialize_params(double *pars)
{
    m_pars = new double[m_Npars];
    std::copy_n(pars,m_Npars, m_pars);
}


/// get parameters
void edgefitting::get_params(double *pars)
{
    std::copy_n(m_pars,m_Npars, pars);
}

/// call the fitting routine, switching betweeen the different lineshape options
void edgefitting::fit(double *x, double *y, int N)
{
    lm_control_struct control = lm_control_double;
    lm_status_struct status;
    control.verbosity = 7;

    printf( "Fitting ...\n" );
    qDebug() << "size of X" << (sizeof(x)/sizeof(*x));
    qDebug() << N;
    switch (myfun){
    case ToFImagingAlgorithms::EdgeTransmissionExponential :
    {
        lmcurve(m_Npars, m_pars, N, x,y, ToFImagingAlgorithms::EdgeFunction::EdgeFunctionTExponential, &control, &status );
        break;
    }
    case ToFImagingAlgorithms::EdgeTransmissionLinear :
    {
        lmcurve(m_Npars, m_pars, N, x,y, ToFImagingAlgorithms::EdgeFunction::EdgeFunctionTLinear, &control, &status );
        break;
    }
    case ToFImagingAlgorithms::EdgeAttenuationExponential :
    {
        lmcurve(m_Npars, m_pars, N, x,y, ToFImagingAlgorithms::EdgeFunction::EdgeFunctionAExponential, &control, &status );
        break;
    }
    case ToFImagingAlgorithms::EdgeAttenuationLinear :
    {
        lmcurve(m_Npars, m_pars, N, x,y, ToFImagingAlgorithms::EdgeFunction::EdgeFunctionALinear, &control, &status );
        break;
    }
    case ToFImagingAlgorithms::EdgeGradientGaussian :
    {
        double *gradient = new double[N];
        kipl::math::num_gradient(y,x,N,gradient);
        lmcurve(m_Npars, m_pars, N, x, gradient, ToFImagingAlgorithms::EdgeFunction::EdgeGradientGaussian, &control, &status);
        delete [] gradient;
        break;
    }
    default :
        throw kipl::base::KiplException("Wrong edge function.",__FILE__,__LINE__);


    }



    printf( "Results:\n" );
    printf( "status after %d function evaluations:\n  %s\n",
            status.nfev, lm_infmsg[status.outcome] );

}

/// compute initial parameters by knowing an estimated Bragg edge position (est_t0)
/// this function is supposed to work only for the complex lineshape formulation (not the gaussian of the gradient)
/// inputs:
/// double *x = x array (ToF or lambda)
/// double *y = actual array
/// int N = length of the array
/// est_t0 estimated value for the edge position (ToF or lambda)
/// This computes the m_pars[3], m_pars[4], m_pars[5], m_pars[6] depending on the lineshape
void edgefitting::compute_initial_params(double *x, double *y, int N, double est_t0)
{
    m_pars[0] = est_t0;
    m_pars[1] = 0.0001; //default?
    m_pars[2] = 0.0015; //default?
}

/// compute initial parameters without knowing an estimated Bragg edge position
/// this function is supposed to work only for the complex lineshape formulation (not the gaussian of the gradient)
/// inputs:
/// double *x = x array (ToF or lambda)
/// double *y = actual array
/// int N = length of the array
/// It first run the gaussian of the gradient to estimate the edge position m_pars[0]
/// Then it computes the m_pars[3], m_pars[4], m_pars[5], m_pars[6] depending on the lineshape
void edgefitting::compute_initial_params(double *x, double *y, int N)
{
    double *gradient = new double[N];
    double *gauss_param = new double[3];
    double *updated_gauss_params = new double[3];
    kipl::math::num_gradient(y,x,N,gradient);

    gauss_param[0] = gauss_param[1] = gauss_param[2] = 1.0;

    ToFImagingAlgorithms::edgefitting myfit(3, ToFImagingAlgorithms::eEdgeFunction::EdgeGradientGaussian);
    myfit.intialize_params(gauss_param);
    myfit.fit(x,y,N);
    myfit.get_params(updated_gauss_params);

    m_pars[0] = updated_gauss_params[0];
    m_pars[1] = 0.0001; //default?
    m_pars[2] = 0.0015; //default?


    delete [] gradient;
    delete [] gauss_param;
    delete [] updated_gauss_params;



}



}
