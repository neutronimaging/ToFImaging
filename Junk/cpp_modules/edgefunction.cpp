#include "edgefunction.h"
#include <strings/miscstring.h>
#include <strings/string2array.h>
#include <base/KiplException.h>
#include <math/mathconstants.h>
#include <math/gradient.h>
#include <cmath>

namespace ToFImagingAlgorithms{

TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::EdgeFunction(int n)
{
    m_Npars=n;
}

TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::~EdgeFunction()
{
    //delete stuff
}

/// \param x The argument
/// \param m_pars The fitting parameter
/// \retval The method returns the lineshape formulation from Santisteban et al. 2001
/// The lineshape is defined by 7 parameters, here:
/// m_pars[0] = t0 is the edge position
/// m_pars[1] = sigma is the Gaussian broadening
/// m_pars[2] = tau is the exponential decay of the trailing edge
/// m_pars[3] = a_{0} first linear parameter for the function after the edge
/// m_pars[4] = b_{0} second linear parameter for the function after the edge
/// m_pars[5] = a_{hkl} first linear parameter for the function before the edge
/// m_pars[6] = b_{hkl} second linear parameter for the function after the edge
///
double TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::EdgeFunctionTExponential(double x, const double *m_pars)
{
    double term3,term4,term5,edge,exp_after,exp_before;
    term3 = std::erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2));
    term4 = std::exp(-1.0*((x-m_pars[0])/m_pars[2])+((m_pars[1]*m_pars[1])/(2.0*m_pars[2]*m_pars[2])));
    term5 = std::erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2)+m_pars[1]/m_pars[2]);
    edge = 0.5*(term3-term4*term5); // myedgefunction
    exp_after = std::exp(-1.0*(m_pars[3]+m_pars[4]*x)); //myexp after
    exp_before = std::exp(-1.0*(m_pars[5]+m_pars[6]*x)); //my exp before
    return exp_after*(exp_before+(1-exp_before)*edge);
}

/// \param x The argument
/// \param m_pars The fitting parameter
/// \retval The method returns the lineshape formulation with linear function modification before and after the edges
/// The lineshape is defined by 7 parameters, here:
/// m_pars[0] = t0 is the edge position
/// m_pars[1] = sigma is the Gaussian broadening
/// m_pars[2] = tau is the exponential decay of the trailing edge
/// m_pars[3] = a_{0} first linear parameter for the function after the edge
/// m_pars[4] = b_{0} second linear parameter for the function after the edge
/// m_pars[5] = a_{hkl} first linear parameter for the function before the edge
/// m_pars[6] = b_{hkl} second linear parameter for the function after the edge
///
double TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::EdgeFunctionTLinear(double x, const double *m_pars)
{
    double term3,term4,term5,edge,line_after,line_before;
    term3 = std::erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2));
    term4 = std::exp(-1.0*((x-m_pars[0])/m_pars[2])+((m_pars[1]*m_pars[1])/(2.0*m_pars[2]*m_pars[2])));
    term5 = std::erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2)+m_pars[1]/m_pars[2]);
    edge = 0.5*(term3-term4*term5); // myedgefunction
    line_after = (m_pars[3]+m_pars[4]*x); //line after
    line_before = (m_pars[5]+m_pars[6]*x); // line before
    return line_after*edge + line_before*(1.0-edge);
}

/// \param x The argument
/// \param m_pars The fitting parameter
/// \retval The method returns the lineshape formulation from Santisteban et al. 2001, with the edge sign modified to accomodate attenuation datasets
/// The lineshape is defined by 7 parameters, here:
/// m_pars[0] = t0 is the edge position
/// m_pars[1] = sigma is the Gaussian broadening
/// m_pars[2] = tau is the exponential decay of the trailing edge
/// m_pars[3] = a_{0} first linear parameter for the function after the edge
/// m_pars[4] = b_{0} second linear parameter for the function after the edge
/// m_pars[5] = a_{hkl} first linear parameter for the function before the edge
/// m_pars[6] = b_{hkl} second linear parameter for the function after the edge
///
double TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::EdgeFunctionAExponential(double x, const double *m_pars)
{
    double term3,term4,term5,edge,exp_after,exp_before;
    term3 = std::erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2));
    term4 = std::exp(-1.0*((x-m_pars[0])/m_pars[2])+((m_pars[1]*m_pars[1])/(2.0*m_pars[2]*m_pars[2])));
    term5 = std::erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2)+m_pars[1]/m_pars[2]);
    edge = 1.0 - 0.5*(term3-term4*term5); // myedgefunction
    exp_after = std::exp(-1.0*(m_pars[3]+m_pars[4]*x)); //myexp after
    exp_before = std::exp(-1.0*(m_pars[5]+m_pars[6]*x)); //my exp before
    return exp_before*(exp_after+(1-exp_after)*edge);
}

/// \param x The argument
/// \param m_pars The fitting parameter
/// \retval The method returns the lineshape formulation with the edge sign modified to accomodate attenuation datasets and linear functions before and after the edge position
/// The lineshape is defined by 7 parameters, here:
/// m_pars[0] = t0 is the edge position
/// m_pars[1] = sigma is the Gaussian broadening
/// m_pars[2] = tau is the exponential decay of the trailing edge
/// m_pars[3] = a_{0} first linear parameter for the function after the edge
/// m_pars[4] = b_{0} second linear parameter for the function after the edge
/// m_pars[5] = a_{hkl} first linear parameter for the function before the edge
/// m_pars[6] = b_{hkl} second linear parameter for the function after the edge
///m_MapPars
double TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::EdgeFunctionALinear(double x, const double *m_pars)
{
    double term3,term4,term5,edge,line_after,line_before;
    term3 = std::erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2));
    term4 = std::exp(-1.0*((x-m_pars[0])/m_pars[2])+((m_pars[1]*m_pars[1])/(2.0*m_pars[2]*m_pars[2])));
    term5 = std::erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2)+m_pars[1]/m_pars[2]);
    edge = 1.0-0.5*(term3-term4*term5); // myedgefunction
    line_after = (m_pars[3]+m_pars[4]*x); // line after
    line_before = (m_pars[5]+m_pars[6]*x); // line before
    return line_before*edge + line_after*(1.0-edge);
}

/// \param x The argument
/// \param m_pars The fitting parameter
/// \retval The method returns the simplified fitting with a Gaussian, the input is expected to be the gradient of the signal
/// The Gaussian is described by the two parameters
/// m_pars[0] = mean, estimating the edge position
/// m_pars[1] = variance, estimating the edge broadening
/// m_pars[2] = amplitude, possibly non useful for the edge fitting
double TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::EdgeGradientGaussian(double x, const double *m_pars)
{
    return m_pars[2]*std::exp(-(x-m_pars[0])*(x-m_pars[0])/(m_pars[1]));
}

/// \param x The argument
/// \param m_pars The fitting parameter
/// \retval The method returns an expontial function with exponent equal to a linear function
/// The exponential is decribed by two parameters as y= exp(-1.0*(m_pars[0]+m_pars[1]*x))
/// m_pars[0] = first parameter
/// m_pars[1] = second parameter
double TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::ExponentialFunction(double x, const double *m_pars)
{
    return std::exp(-1.0*(m_pars[0]+m_pars[1]*x));
}

/// \param x The argument
/// \param m_pars The fitting parameter
/// \retval The method returns an expontial function with exponent equal to a linear function
/// The exponential is decribed by two parameters as y= exp(-1.0*(m_pars[0]+m_pars[1]*x))
/// m_pars[0] = first parameter
/// m_pars[1] = second parameter
/// m_pars[2] = third parameter
/// m_pars[3] = fourth parameter
double TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::CombinedExponentialFunction(double x, const double *m_pars)
{
    double *m_pars_after = new double[2];
    m_pars_after[0] = m_pars[0];
    m_pars_after[1] = m_pars[1];

    double *m_pars_before = new double [2];
    m_pars_before[0] = m_pars[2];
    m_pars_before[1] = m_pars[3];

    return ExponentialFunction(x,m_pars_after)*ExponentialFunction(x,m_pars_before);
}

TOF_IMAGINGALGORITHMSHARED_EXPORT void string2enum(std::string &str, ToFImagingAlgorithms::eEdgeFunction &e)
{
    std::string lowstr=kipl::strings::toLower(str);

    std::map<std::string,ToFImagingAlgorithms::eEdgeFunction> convmap;

    convmap["edgetlinear"] = ToFImagingAlgorithms::eEdgeFunction::EdgeTransmissionLinear;
    convmap["edgetexponential"] = ToFImagingAlgorithms::eEdgeFunction::EdgeTransmissionExponential;
    convmap["edgealinear"] = ToFImagingAlgorithms::eEdgeFunction::EdgeAttenuationLinear;
    convmap["edgeaexponential"] = ToFImagingAlgorithms::eEdgeFunction::EdgeAttenuationExponential;
    convmap["edgegradgaussian"] = ToFImagingAlgorithms::eEdgeFunction::EdgeGradientGaussian;

    auto it=convmap.find(lowstr);

    if (it==convmap.end())
        throw kipl::base::KiplException("Profile function does not exist",__FILE__,__LINE__);

    e=it->second;

}

TOF_IMAGINGALGORITHMSHARED_EXPORT std::string enum2string(ToFImagingAlgorithms::eEdgeFunction  e)
{
    std::string str;
    switch(e){
        case ToFImagingAlgorithms::EdgeTransmissionLinear:
            str = "EdgeTLinear";
            break;
        case ToFImagingAlgorithms::EdgeTransmissionExponential:
            str =  "EdgeTExponential";
            break;
        case ToFImagingAlgorithms::EdgeAttenuationLinear:
            str =  "EdgeALinear";
            break;
        case ToFImagingAlgorithms::EdgeAttenuationExponential:
            str =  "EdgeAExponential";
            break;
        case ToFImagingAlgorithms::EdgeGradientGaussian:
            str =  "EdgeGradGaussian";
            break;
    default:
        str =  "EdgeTExponential"; // or throw Exception

    }
    return str;
}
TOF_IMAGINGALGORITHMSHARED_EXPORT std::ostream  &  operator<<(std::ostream &s,ToFImagingAlgorithms::eEdgeFunction e)
{
    s<<enum2string(e);
    return s;
}


}
