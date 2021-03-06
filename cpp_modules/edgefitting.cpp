#include "edgefitting.h"
#include <findclosest.h>
#include <base/KiplException.h>
#include <lmcurve.h>
#include <QDebug>
#include <math/gradient.h>
#include <math/linfit.h>
#include <iterator>
#include <filters/savitzkygolayfilter.h>

namespace ToFImagingAlgorithms {

/// class constructor, inputs:
/// \param n number of parameters to be fitted, 7 for the complete lineshape, 3 for the simplifies (gaussian of the signal gradient)
/// \param ef lineshape type
TOF_IMAGINGALGORITHMSHARED_EXPORT edgefitting::edgefitting(int n, ToFImagingAlgorithms::eEdgeFunction ef):
    m_Npars(n),
    myfun(ef),
    blinear(true),
    bsmooth(false),
    m_sgWin(5),
    m_sgPoly(2)
{

    m_pars.assign(m_Npars,1.0); // initialize container


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
TOF_IMAGINGALGORITHMSHARED_EXPORT edgefitting::~edgefitting()
{

}


/// compute initial parameters without knowing an estimated Bragg edge position
/// this function is supposed to work only for the complex lineshape formulation (not the gaussian of the gradient)
/// inputs:
/// \param x x-array for the edge signal (ToF or lambda)
/// \param y edge signal
/// \param N length of the array
/// It first run the gaussian of the gradient to estimate the edge position m_pars[0]
/// Then it computes the m_pars[3], m_pars[4], m_pars[5], m_pars[6] depending on the lineshape
void TOF_IMAGINGALGORITHMSHARED_EXPORT edgefitting::compute_initial_params(std::vector<double> &x, std::vector<double> &y, int N)
{

    unsigned int est_pos, size_1, size_2, buffer;
    std::vector <double> gauss_param(3), updated_gauss_params(3);
    ToFImagingAlgorithms::edgefitting myfit(3, ToFImagingAlgorithms::eEdgeFunction::EdgeGradientGaussian);
    myfit.fit(x,y,N);
    myfit.get_params(updated_gauss_params);

    m_pars[0] = updated_gauss_params[0];
    m_pars[1] = 0.0001; //default?
    m_pars[2] = 0.0015; //default?

//    std::cout << updated_gauss_params[0] << std::endl;

    buffer = static_cast<unsigned int>(0.1*N);
    est_pos = ToFImagingAlgorithms::findClosest(&(x[0]),N, m_pars[0]);


    size_1 = est_pos-buffer;
    size_2 = N-(est_pos+buffer);

    std::vector<double> x1(size_1), y1(size_1), x2(size_2), y2(size_2), logy1(size_1), logy2(size_2);

    std::copy_n(x.begin(), size_1, x1.begin());
    std::copy_n(y.begin(), size_1, y1.begin());
    std::copy_n(x.begin()+(est_pos+buffer), size_2, x2.begin());
    std::copy_n(y.begin()+(est_pos+buffer), size_2, y2.begin());

    double lin_par_before[2];
    double lin_par_after[2];

    if (blinear){
        kipl::math::LinearLSFit(&(x2[0]),&(y2[0]),size_2, lin_par_after, lin_par_after+1, nullptr);
        kipl::math::LinearLSFit(&(x1[0]),&(y1[0]),size_1, lin_par_before, lin_par_before+1, nullptr);
    }
    else {

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




}

/// compute initial parameters by knowing an estimated Bragg edge position (est_t0)
/// this function is supposed to work only for the complex lineshape formulation (not the gaussian of the gradient)
/// inputs:
/// \param x x-array for the edge signal (ToF or lambda)
/// \param y edge signal
/// \param N length of the array
/// \param est_t0 estimated value for the edge position (ToF or lambda)
/// This computes the m_pars[3], m_pars[4], m_pars[5], m_pars[6] depending on the lineshape
void TOF_IMAGINGALGORITHMSHARED_EXPORT edgefitting::compute_initial_params(std::vector<double> &x, std::vector<double> &y, int N, double est_t0)
{
    unsigned int est_pos, size_1, size_2, buffer;
    m_pars[0] = est_t0;
    m_pars[1] = 0.0001; //default?
    m_pars[2] = 0.0015; //default?

    buffer = static_cast<unsigned int>(0.1*N); // this is possibly to optimize
    est_pos = ToFImagingAlgorithms::findClosest(&(x[0]),N, est_t0);
    size_1 = est_pos-buffer;
    size_2 = N-(est_pos+buffer);

    std::vector<double> x1(size_1),x2(size_2),y1(size_1),y2(size_2),logy1(size_1),logy2(size_2);

    // divide the signal in two


    std::copy_n(x.begin(), size_1, x1.begin());
    std::copy_n(y.begin(), size_1, y1.begin());
    std::copy_n(x.begin()+(est_pos+buffer), size_2, x2.begin());
    std::copy_n(y.begin()+(est_pos+buffer), size_2, y2.begin());


    double lin_par_before[2];
    double lin_par_after[2];

    if (blinear){
        kipl::math::LinearLSFit(&(x2[0]),&(y2[0]),size_2, lin_par_after, lin_par_after+1, nullptr);
        kipl::math::LinearLSFit(&(x1[0]),&(y1[0]),size_1, lin_par_before, lin_par_before+1, nullptr);
    }
    else {

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

}



/// call the fitting routine, switching betweeen the different lineshape options
/// \param x x-axis for the edge signal
/// \param y edge signal
void TOF_IMAGINGALGORITHMSHARED_EXPORT edgefitting::fit(std::vector<double> &x, std::vector<double> &y, int N)
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
        std::vector<double> gradient(N),smoothed(N);
        kipl::math::num_gradient(&(y[0]),&(x[0]),N,&(gradient[0]));
        lmcurve(m_Npars, &(m_pars[0]), N, &(x[0]), &(gradient[0]), ToFImagingAlgorithms::EdgeFunction::EdgeGradientGaussian, &control, &status);
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

/// function that the returns the fitting parameters
/// \param pars the parameters where the par values are written
void TOF_IMAGINGALGORITHMSHARED_EXPORT edgefitting::get_params(std::vector<double> &pars)
{
    pars.assign(m_pars.begin(), m_pars.end());
}

/// initialize parameters by copying the inputs in the m_pars
/// \param pars fitting parameters
void TOF_IMAGINGALGORITHMSHARED_EXPORT edgefitting::initialize_params(std::vector<double> &pars)
{
    m_pars.assign(pars.begin(), pars.end());
}


/// apply sav-gol filter to smooth the signal
/// inputs:
/// \param x x axis for the edge (currently not used)
/// \param y signal to be smoothed (edge)
void TOF_IMAGINGALGORITHMSHARED_EXPORT edgefitting::smooth(std::vector<double> &x, std::vector<double> &y)
{
    kipl::filters::SavitzkyGolayFilter<double> sg;
    y = sg(y,m_sgWin,m_sgPoly); // input is same as output. is it ok?

}






}
