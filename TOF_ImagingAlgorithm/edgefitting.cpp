#include "edgefitting.h"
#include <findclosest.h>
#include <base/KiplException.h>
#include <lmcurve.h>
#include <QDebug>
#include <math/gradient.h>
#include <math/linfit.h>
#include <iterator>

namespace ToFImagingAlgorithms {

/// class constructor, inputs:
/// n = number of parameters to be fitted, 7 for the complete lineshape, 3 for the simplifies (gaussian of the signal gradient)
/// ef = lineshape type
edgefitting::edgefitting(int n, ToFImagingAlgorithms::eEdgeFunction ef)
{
    m_Npars = n;
    myfun = ef;
    blinear = true;

//    EdgeTransmissionExponential,
//    EdgeTransmissionLinear,
//    EdgeAttenuationExponential,
//    EdgeAttenuationLinear,
//    EdgeGradientGaussian,

    switch (myfun) {
        case ToFImagingAlgorithms::EdgeTransmissionLinear :
        {
            blinear = true;
            break;
        }

        case ToFImagingAlgorithms::EdgeAttenuationLinear :
        {
            blinear = true;
            break;
        }

        case  ToFImagingAlgorithms::EdgeTransmissionExponential :
        {
            blinear = false;
            break;
        }

        case ToFImagingAlgorithms::EdgeAttenuationExponential :
        {
            blinear = false;
            break;
        }

        case ToFImagingAlgorithms::EdgeGradientGaussian :
        {
//            blinear = false; // this is not relavant in this case
            break;
        }

    }

}

/// class deconstructor
edgefitting::~edgefitting()
{

}


/// initialize parameters by copying the inputs in the m_pars
void edgefitting::intialize_params(std::vector<double> &pars)
{   
    m_pars.assign(pars.begin(), pars.end());
}

/// get parameters
void edgefitting::get_params(std::vector<double> &pars)
{
    pars.assign(m_pars.begin(), m_pars.end());
}

/// call the fitting routine, switching betweeen the different lineshape options
void edgefitting::fit(std::vector<double> &x, std::vector<double> &y, int N)
{
    lm_control_struct control = lm_control_double;
    lm_status_struct status;
    control.verbosity = 7;

    printf( "Fitting ...\n" );
    qDebug() << "size of X" << x.size();
    qDebug() << N;
    switch (myfun){
    case ToFImagingAlgorithms::EdgeTransmissionExponential :
    {
        lmcurve(m_Npars, &(m_pars[0]), N, &(x[0]),&(y[0]), ToFImagingAlgorithms::EdgeFunction::EdgeFunctionTExponential, &control, &status );
        break;
    }
    case ToFImagingAlgorithms::EdgeTransmissionLinear :
    {
        lmcurve(m_Npars, &(m_pars[0]), N, &(x[0]),&(y[0]), ToFImagingAlgorithms::EdgeFunction::EdgeFunctionTLinear, &control, &status );
        break;
    }
    case ToFImagingAlgorithms::EdgeAttenuationExponential :
    {
        lmcurve(m_Npars, &(m_pars[0]), N, &(x[0]),&(y[0]), ToFImagingAlgorithms::EdgeFunction::EdgeFunctionAExponential, &control, &status );
        break;
    }
    case ToFImagingAlgorithms::EdgeAttenuationLinear :
    {
        lmcurve(m_Npars, &(m_pars[0]), N, &(x[0]),&(y[0]), ToFImagingAlgorithms::EdgeFunction::EdgeFunctionALinear, &control, &status );
        break;
    }
    case ToFImagingAlgorithms::EdgeGradientGaussian :
    {
        double *gradient = new double[N];
        kipl::math::num_gradient(&(y[0]),&(x[0]),N,gradient);
        lmcurve(m_Npars, &(m_pars[0]), N, &(x[0]), gradient, ToFImagingAlgorithms::EdgeFunction::EdgeGradientGaussian, &control, &status);
        delete [] gradient;
        break;
    }
    default :
    {
        throw kipl::base::KiplException("Wrong edge function.",__FILE__,__LINE__);
    }


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
void edgefitting::compute_initial_params(std::vector<double> &x, std::vector<double> &y, int N, double est_t0)
{
    unsigned int est_pos, size_1, size_2, buffer;

//    double *x1, *x2, *y1, *y2, *logy1, *logy2;
    m_pars[0] = est_t0;
    m_pars[1] = 0.0001; //default?
    m_pars[2] = 0.0015; //default?

    buffer = static_cast<unsigned int>(0.1*N); // this is possibly to optimize
    est_pos = ToFImagingAlgorithms::findClosest(&(x[0]),N, est_t0);
    size_1 = est_pos-buffer;
    size_2 = N-(est_pos+buffer);

    std::vector<double> x1(size_1),x2(size_2),y1(size_1),y2(size_2),logy1(size_1),logy2(size_2);

    // divide the signal in two

//    x1 = new double[size_1];
//    y1 = new double[size_1];
//    x2 = new double[size_2];
//    y2 = new double[size_2];

    std::copy_n(x.begin(), size_1, std::back_inserter(x1));
    std::copy_n(y.begin(), size_1, std::back_inserter(y1));
    std::copy_n(x.at(est_pos+buffer), size_2, std::back_inserter(x2));
    std::copy_n(y.at(est_pos+buffer), size_2, std::back_inserter(y2));


    double lin_par_before[2];
    double lin_par_after[2];

    if (blinear){
        kipl::math::LinearLSFit(&(x2[0]),&(y2[0]),size_2, lin_par_after, lin_par_after+1, nullptr);
        kipl::math::LinearLSFit(&(x1[0]),&(y1[0]),size_1, lin_par_before, lin_par_before+1, nullptr);
    }
    else {
//        logy1 = new double[size_1];
//        logy2 = new double[size_2];

        //compute the log of the ys
        for (unsigned int i=0; i<size_2; ++i){
            logy2[i] = -1.0*std::log(y2[i]);
        }


        kipl::math::LinearLSFit(&(x2[0]), &(logy2[0]), size_2, lin_par_after, lin_par_after+1, nullptr);

        for (unsigned int i=0; i<size_1; ++i)
        {
            logy1[i] = -1.0*std::log(y1[i])-lin_par_after[0]-lin_par_after[1]*x[0];
        }

        kipl::math::LinearLSFit(&(x1[0]), &(logy1[0]), size_1, lin_par_before, lin_par_before+1, nullptr);


    }



    m_pars[3] = lin_par_after[0];
    m_pars[4] = lin_par_after[1];
    m_pars[5] = lin_par_before[0];
    m_pars[6] = lin_par_before[1];

//    delete [] x1;
//    delete [] x2;
//    delete [] y1;
//    delete [] y2;


}

/// compute initial parameters without knowing an estimated Bragg edge position
/// this function is supposed to work only for the complex lineshape formulation (not the gaussian of the gradient)
/// inputs:
/// double *x = x array (ToF or lambda)
/// double *y = actual array
/// int N = length of the array
/// It first run the gaussian of the gradient to estimate the edge position m_pars[0]
/// Then it computes the m_pars[3], m_pars[4], m_pars[5], m_pars[6] depending on the lineshape
void edgefitting::compute_initial_params(std::vector<double> &x, std::vector<double> &y, int N)
{

     unsigned int est_pos, size_1, size_2, buffer;
//    double *x1, *x2, *y1, *y2, *logy1, *logy2;


//    double *gradient = new double[N];
//    double *gauss_param = new double[3];
//    double *updated_gauss_params = new double[3];
    std::vector <double> gradient, gauss_param, updated_gauss_params;
    kipl::math::num_gradient(&(y[0]),&(x[0]),N,&(gradient[0]));

    gauss_param[0] = gauss_param[1] = gauss_param[2] = 1.0;

    ToFImagingAlgorithms::edgefitting myfit(3, ToFImagingAlgorithms::eEdgeFunction::EdgeGradientGaussian);
    myfit.intialize_params(gauss_param);
    myfit.fit(x,y,N);
    myfit.get_params(updated_gauss_params);

    m_pars[0] = updated_gauss_params[0];
    m_pars[1] = 0.0001; //default?
    m_pars[2] = 0.0015; //default?

    std::cout << m_pars[0] << std::endl;

    buffer = static_cast<unsigned int>(0.1*N);
    est_pos = ToFImagingAlgorithms::findClosest(&(x[0]),N, m_pars[0]);


    size_1 = est_pos-buffer;
    size_2 = N-(est_pos+buffer);

    std::vector<double> x1(size_1), y1(size_1), x2(size_2), y2(size_2), logy1(size_1), logy2(size_2);

//    x1 = new double[size_1];
//    y1 = new double[size_1];
//    x2 = new double[size_2];
//    y2 = new double[size_2];

    std::copy_n(x.begin(), size_1, std::back_inserter(x1));
    std::copy_n(y.begin(), size_1, std::back_inserter(y1));
    std::copy_n(x.at(est_pos+buffer), size_2, std::back_inserter(x2));
    std::copy_n(y.at(est_pos+buffer), size_2, std::back_inserter(y2));

    double lin_par_before[2];
    double lin_par_after[2];

    if (blinear){
        kipl::math::LinearLSFit(&(x2[0]),&(y2[0]),size_2, lin_par_after, lin_par_after+1, nullptr);
        kipl::math::LinearLSFit(&(x1[0]),&(y1[0]),size_1, lin_par_before, lin_par_before+1, nullptr);
    }
    else {
//        logy1 = new double[size_1];
//        logy2 = new double[size_2];

        //compute the log of the ys
        for (unsigned int i=0; i<size_2; ++i){
            logy2[i] = -1.0*std::log(y2[i]);
        }

        kipl::math::LinearLSFit(&(x2[0]), &(logy2[0]), size_2, lin_par_after, lin_par_after+1, nullptr);

        for (unsigned int i=0; i<size_1; ++i)
        {
            logy1[i] = -1.0*std::log(y1[i])-lin_par_after[0]-lin_par_after[1]*x[0];
        }

        kipl::math::LinearLSFit(&(x1[0]), &(logy1[0]), size_1, lin_par_before, lin_par_before+1, nullptr);


    }


    m_pars[3] = lin_par_after[0];
    m_pars[4] = lin_par_after[1];
    m_pars[5] = lin_par_before[0];
    m_pars[6] = lin_par_before[1];

//    delete [] x1;
//    delete [] x2;
//    delete [] y1;
//    delete [] y2;


//    delete [] gradient;
//    delete [] gauss_param;
//    delete [] updated_gauss_params;



}



}
