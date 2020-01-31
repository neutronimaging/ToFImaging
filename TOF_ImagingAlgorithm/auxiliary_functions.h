#ifndef AUXILIARY_FUNCTIONS_H
#define AUXILIARY_FUNCTIONS_H

namespace ToFImagingAlgorithms {

///\brief linear function
double LinearFunction(double x, const double *m_pars); // possibly not used
///\brief exponential function
double ExponentialFunction(double x, const double *m_pars);
///\brief product of exponentials function
double CombinedExponentialFunction(double x, const double *m_pars);

}

#endif // AUXILIARY_FUNCTIONS_H
