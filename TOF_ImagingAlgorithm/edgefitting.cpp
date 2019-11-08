#include "edgefitting.h"
#include <lmcurve.h>
#include <QDebug>
#include <math/gradient.h>
#include <base/KiplException.h>

edgefitting::edgefitting(int n, BraggEdge::eEdgeFunction ef)
{
    m_Npars = n;
    myfun = ef;
}

edgefitting::~edgefitting()
{

}

void edgefitting::intialize_params(double *pars)
{
    m_pars = new double[m_Npars];
    std::copy_n(pars,m_Npars, m_pars);
}

void edgefitting::get_params(double *pars)
{
    std::copy_n(m_pars,m_Npars, pars);
}

void edgefitting::fit(double *x, double *y, int N)
{
    lm_control_struct control = lm_control_double;
    lm_status_struct status;
    control.verbosity = 7;

    printf( "Fitting ...\n" );
    qDebug() << "size of X" << (sizeof(x)/sizeof(*x));
    qDebug() << N;
    switch (myfun){
    case BraggEdge::EdgeTransmissionExponential :
    {
        lmcurve(m_Npars, m_pars, N, x,y, BraggEdge::EdgeFunction::EdgeFunctionTExponential, &control, &status );
        break;
    }
    case BraggEdge::EdgeTransmissionLinear :
    {
        lmcurve(m_Npars, m_pars, N, x,y, BraggEdge::EdgeFunction::EdgeFunctionTLinear, &control, &status );
        break;
    }
    case BraggEdge::EdgeAttenuationExponential :
    {
        lmcurve(m_Npars, m_pars, N, x,y, BraggEdge::EdgeFunction::EdgeFunctionAExponential, &control, &status );
        break;
    }
    case BraggEdge::EdgeAttenuationLinear :
    {
        lmcurve(m_Npars, m_pars, N, x,y, BraggEdge::EdgeFunction::EdgeFunctionALinear, &control, &status );
        break;
    }
    case BraggEdge::EdgeGradientGaussian :
    {
        double *gradient = new double[N];
        kipl::math::num_gradient(y,x,N,gradient);
        lmcurve(m_Npars, m_pars, N, x, gradient, BraggEdge::EdgeFunction::EdgeGradientGaussian, &control, &status);
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
