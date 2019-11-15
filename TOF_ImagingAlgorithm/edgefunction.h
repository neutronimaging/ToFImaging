#ifndef EDGEFUNCTION_H
#define EDGEFUNCTION_H

#include "tof_imagingalgorithm_global.h"
#include <string>
#include <map>


namespace ToFImagingAlgorithms{
/// enum of different edge functions available
enum eEdgeFunction{
    EdgeTransmissionExponential,
    EdgeTransmissionLinear,
    EdgeAttenuationExponential,
    EdgeAttenuationLinear,
    EdgeGradientGaussian,
};



class TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction
{
public:
    EdgeFunction(int n=1);
    ~EdgeFunction();
    ///\brief Implements edge function in transmission with exponentials before and after the edge
    static double EdgeFunctionTExponential(double x, const double *m_pars);
    ///\brief Implements edge function in transmission with lines before and after the edge
    static double EdgeFunctionTLinear(double x, const double *m_pars);
    ///\brief Implements edge function in attenuation with lines before and after the edge
    static double EdgeFunctionAExponential(double x, const double *m_pars);
    ///\brief Implements edge function in attenuation with lines before and after the edge
    static double EdgeFunctionALinear(double x, const double *m_pars);
    ///\brief Implements the Gaussian shape for simple fitting of the edge position (Gaussian fitting of the first derivative of the signal)
    static double EdgeGradientGaussian(double x, const double *m_pars);




protected:
//    /// \brief Parameter array
//    double *m_pars;
//    /// \brief Parameter lock array
//    bool *m_lock;
    /// \brief The number of parameters
    int m_Npars;
//    /// \brief The number of parameters to be fitted
//    int m_pars2fit;
};

}

void  TOF_IMAGINGALGORITHMSHARED_EXPORT string2enum(std::string &str, ToFImagingAlgorithms::eEdgeFunction &e);
TOF_IMAGINGALGORITHMSHARED_EXPORT std::string enum2string(ToFImagingAlgorithms::eEdgeFunction  e);
TOF_IMAGINGALGORITHMSHARED_EXPORT std::ostream  &  operator<<(std::ostream &s,ToFImagingAlgorithms::eEdgeFunction e);

#endif // EDGEFUNCTION_H
