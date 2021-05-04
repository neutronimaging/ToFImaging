#ifndef AUXILIARY_FUNCTIONS_H
#define AUXILIARY_FUNCTIONS_H

#include "tof_imagingalgorithm_global.h"

namespace ToFImagingAlgorithms {

///\brief linear function
double TOF_IMAGINGALGORITHMSHARED_EXPORT LinearFunction(double x, const double *m_pars); // possibly not used
///\brief exponential function
double TOF_IMAGINGALGORITHMSHARED_EXPORT ExponentialFunction(double x, const double *m_pars);
///\brief product of exponentials function
double TOF_IMAGINGALGORITHMSHARED_EXPORT CombinedExponentialFunction(double x, const double *m_pars);

}

#endif // AUXILIARY_FUNCTIONS_H
