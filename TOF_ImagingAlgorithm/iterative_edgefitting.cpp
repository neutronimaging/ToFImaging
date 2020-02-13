#include "iterative_edgefitting.h"

namespace ToFImagingAlgorithms {

/// class constructor
TOF_IMAGINGALGORITHMSHARED_EXPORT iterative_edgefitting::iterative_edgefitting(int n, ToFImagingAlgorithms::eEdgeFunction ef):
    m_Npars(n),
    myfun(ef),
    edgefitting (m_Npars, myfun)
{

}

/// class destructor
TOF_IMAGINGALGORITHMSHARED_EXPORT iterative_edgefitting::~iterative_edgefitting()
{

}

}
