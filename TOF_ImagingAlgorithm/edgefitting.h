#ifndef EDGEFITTING_H
#define EDGEFITTING_H

#include <edgefunction.h>
namespace  ToFImagingAlgorithms

{
class edgefitting
{
public:
    edgefitting(int n, ToFImagingAlgorithms::eEdgeFunction ef);
    void intialize_params(double *pars);
    void get_params(double *pars);
    void fit(double *x, double *y, int N);
    ~edgefitting();


private:
//    BraggEdge::EdgeFunction edge_lineshape(); /// lineshape to be fitted // does not have to be a class parameter
//    double *x; /// x  axis of the line to be fitted // don't need them either probably
//    double *y; /// y  axis of the line to be fitted
    double *m_pars; /// parameter arrays
    int m_Npars; /// number of parameters
    ToFImagingAlgorithms::eEdgeFunction myfun; /// type of function to be fitted
};


}


#endif // EDGEFITTING_H
