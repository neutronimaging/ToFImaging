#include "auxiliary_functions.h"
#include <math.h>

namespace ToFImagingAlgorithms {

double LinearFunction(double x, const double *m_pars)
{
    return (m_pars[0]+m_pars[1]*x);
}

double ExponentialFunction(double x, const double *m_pars)
{
    return exp(-1.0*(m_pars[0]+m_pars[1]*x));
}

double CombinedExponentialFunction(double x, const double *m_pars)
{
    return exp(-1.0*(m_pars[0]+m_pars[1]*x))*exp(-1.0*(m_pars[3]+m_pars[4]*x));
}

}
