#include "edgefitting.h"
#include <lmcurve.h>
#include <QDebug>

edgefitting::edgefitting(int n)
{
    m_Npars = n;
}

edgefitting::~edgefitting()
{

}

void edgefitting::intialize_params(double *pars)
{
    m_pars = new double[m_Npars];
    m_pars = pars;
}

double *edgefitting::get_params()
{
    return m_pars;
}

void edgefitting::fit(double *x, double *y, int N)
{
    lm_control_struct control = lm_control_double;
    lm_status_struct status;
    control.verbosity = 7;

    printf( "Fitting ...\n" );
    qDebug() << "size of X" << (sizeof(x)/sizeof(*x));
    qDebug() << N;
    lmcurve(m_Npars, m_pars, N, x,y, BraggEdge::EdgeFunction::EdgeFunctionTExponential, &control, &status );

    printf( "Results:\n" );
    printf( "status after %d function evaluations:\n  %s\n",
            status.nfev, lm_infmsg[status.outcome] );

}
