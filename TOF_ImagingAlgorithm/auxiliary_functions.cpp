#include "auxiliary_functions.h"
#include <math.h>

namespace ToFImagingAlgorithms {

double TOF_IMAGINGALGORITHMSHARED_EXPORT LinearFunction(double x, const double *m_pars)
{
    return (m_pars[0]+m_pars[1]*x);
}

double TOF_IMAGINGALGORITHMSHARED_EXPORT ExponentialFunction(double x, const double *m_pars)
{
    return exp(-1.0*(m_pars[0]+m_pars[1]*x));
}

double TOF_IMAGINGALGORITHMSHARED_EXPORT CombinedExponentialFunction(double x, const double *m_pars)
{
    return exp(-1.0*(m_pars[0]+m_pars[1]*x))*exp(-1.0*(m_pars[3]+m_pars[4]*x));
}

}
