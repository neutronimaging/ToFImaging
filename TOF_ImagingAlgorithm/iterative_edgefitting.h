#ifndef ITERATIVE_EDGEFITTING_H
#define ITERATIVE_EDGEFITTING_H

#include "tof_imagingalgorithm_global.h"
#include <edgefunction.h>
#include <edgefitting.h>


namespace  ToFImagingAlgorithms

{
class iterative_edgefitting : public edgefitting
{
public:
    iterative_edgefitting(int n, ToFImagingAlgorithms::eEdgeFunction ef);
    ~iterative_edgefitting();

private:
   std::map<std::string, double> m_MapPars; /// \brief map for the parmeters
   int m_Npars; ///\param m_Npars number of parameters
   ToFImagingAlgorithms::eEdgeFunction myfun; ///\param myfum type of function to be fitted
};

}
#endif // ITERATIVE_EDGEFITTING_H
