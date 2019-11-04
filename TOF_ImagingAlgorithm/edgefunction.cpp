#include "edgefunction.h"
#include <strings/miscstring.h>
#include <strings/string2array.h>
#include <base/KiplException.h>
#include <math/mathconstants.h>
#include <math/gradient.h>
#include <math.h>

namespace BraggEdge{

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
    term3 = erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2));
    term4 = exp(-1.0*((x-m_pars[0])/m_pars[2])+((m_pars[1]*m_pars[1])/(2.0*m_pars[2]*m_pars[2])));
    term5 = erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2)+m_pars[1]/m_pars[2]);
    edge = 0.5*(term3-term4*term5); // myedgefunction
    exp_after = exp(-1.0*(m_pars[3]+m_pars[4]*x)); //myexp after
    exp_before = exp(-1.0*(m_pars[5]+m_pars[6]*x)); //my exp before
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
    term3 = erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2));
    term4 = exp(-1.0*((x-m_pars[0])/m_pars[2])+((m_pars[1]*m_pars[1])/(2.0*m_pars[2]*m_pars[2])));
    term5 = erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2)+m_pars[1]/m_pars[2]);
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
    term3 = erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2));
    term4 = exp(-1.0*((x-m_pars[0])/m_pars[2])+((m_pars[1]*m_pars[1])/(2.0*m_pars[2]*m_pars[2])));
    term5 = erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2)+m_pars[1]/m_pars[2]);
    edge = 1.0 - 0.5*(term3-term4*term5); // myedgefunction
    exp_after = exp(-1.0*(m_pars[3]+m_pars[4]*x)); //myexp after
    exp_before = exp(-1.0*(m_pars[5]+m_pars[6]*x)); //my exp before
    return exp_after*(exp_before+(1-exp_before)*edge);
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
///
double TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::EdgeFunctionALinear(double x, const double *m_pars)
{
    double term3,term4,term5,edge,line_after,line_before;
    term3 = erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2));
    term4 = exp(-1.0*((x-m_pars[0])/m_pars[2])+((m_pars[1]*m_pars[1])/(2.0*m_pars[2]*m_pars[2])));
    term5 = erfc(-1.0*(x-m_pars[0])/(m_pars[1]*dsqrt2)+m_pars[1]/m_pars[2]);
    edge = 1.0-0.5*(term3-term4*term5); // myedgefunction
    line_after = (m_pars[3]+m_pars[4]*x); // line after
    line_before = (m_pars[5]+m_pars[6]*x); // line before
    return line_after*edge + line_before*(1.0-edge);
}

/// \param x The argument
/// \param m_pars The fitting parameter
/// \retval The method returns the simplified fitting with a Gaussian, the input is expected to be the gradient of the signal
/// The Gaussian is described by the two parameters
/// m_pars[0] = mean, estimating the edge position
/// m_pars[1] = standard deviation, estimating the edge broadening
double TOF_IMAGINGALGORITHMSHARED_EXPORT EdgeFunction::EdgeGradientGaussian(double x, const double *m_pars)
{
    return exp(-(x-m_pars[0])*(x-m_pars[0])/(2.0*m_pars[1]*m_pars[1]));
}

TOF_IMAGINGALGORITHMSHARED_EXPORT void string2enum(std::string &str, BraggEdge::eEdgeFunction &e)
{
    std::string lowstr=kipl::strings::toLower(str);

    std::map<std::string,BraggEdge::eEdgeFunction> convmap;

    convmap["edgetlinear"] = BraggEdge::eEdgeFunction::EdgeTransmissionLinear;
    convmap["edgetexponential"] = BraggEdge::eEdgeFunction::EdgeTransmissionExponential;
    convmap["edgealinear"] = BraggEdge::eEdgeFunction::EdgeAttenuationLinear;
    convmap["edgeaexponential"] = BraggEdge::eEdgeFunction::EdgeAttenuationExponential;
    convmap["edgegradgaussian"] = BraggEdge::eEdgeFunction::EdgeGradientGaussian;

    auto it=convmap.find(lowstr);

    if (it==convmap.end())
        throw kipl::base::KiplException("Profile function does not exist",__FILE__,__LINE__);

    e=it->second;

}

TOF_IMAGINGALGORITHMSHARED_EXPORT std::string enum2string(BraggEdge::eEdgeFunction  e)
{
    std::string str;
    switch(e){
        case BraggEdge::EdgeTransmissionLinear:
            str = "EdgeTLinear";
            break;
        case BraggEdge::EdgeTransmissionExponential:
            str =  "EdgeTExponential";
            break;
        case BraggEdge::EdgeAttenuationLinear:
            str =  "EdgeALinear";
            break;
        case BraggEdge::EdgeAttenuationExponential:
            str =  "EdgeAExponential";
            break;
        case BraggEdge::EdgeGradientGaussian:
            str =  "EdgeGradGaussian";
            break;
    default:
        str =  "EdgeTExponential"; // or throw Exception

    }
    return str;
}
TOF_IMAGINGALGORITHMSHARED_EXPORT std::ostream  &  operator<<(std::ostream &s,BraggEdge::eEdgeFunction e)
{
    s<<enum2string(e);
    return s;
}


}