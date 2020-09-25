#include "iterative_edgefitting.h"

#include <ImagingException.h>

namespace ToFImagingAlgorithms {

/// class constructor
TOF_IMAGINGALGORITHMSHARED_EXPORT iterative_edgefitting::iterative_edgefitting(int n, ToFImagingAlgorithms::eEdgeFunction ef):
    edgefitting (n, ef),
    m_Npars(n),
    myfun(ef)

{
    if (m_Npars==7) // iterative fitting is only relevant for the 7 parameters edge lineshape
    {
        m_MapPars["t0"] = 1.0;
        m_MapPars["sigma"] = 1.0;
        m_MapPars["tau"] = 1.0;
        m_MapPars["a_0"] = 1.0;
        m_MapPars["b_0"] = 1.0;
        m_MapPars["a_hkl"] = 1.0;
        m_MapPars["b_hkl"] = 1.0;
    }
    else {
        throw ImagingException("Wrong number of parameters", __FILE__, __LINE__);
    }

}

/// class destructor
TOF_IMAGINGALGORITHMSHARED_EXPORT iterative_edgefitting::~iterative_edgefitting()
{

}

void TOF_IMAGINGALGORITHMSHARED_EXPORT iterative_edgefitting::fit()
{
    //fit
}

void TOF_IMAGINGALGORITHMSHARED_EXPORT iterative_edgefitting::initialize_map(double est_pos)
{
    m_MapPars["t0"] = est_pos;
}

void TOF_IMAGINGALGORITHMSHARED_EXPORT iterative_edgefitting::initialize_map(std::vector<double> &param)
{
    m_MapPars["t0"] = param[0];
    m_MapPars["sigma"] = param[1];
    m_MapPars["tau"] = param[2];
    m_MapPars["a_0"] = param[3];
    m_MapPars["b_0"] = param[4];
    m_MapPars["a_hkl"] = param[5];
    m_MapPars["b_hkl"] = param[6];
}

void TOF_IMAGINGALGORITHMSHARED_EXPORT iterative_edgefitting::set_fixedparam(std::vector<std::string> &fixed_param)
{
    // set the fixed params
}

}
