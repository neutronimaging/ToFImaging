#ifndef EDGEFITTING_H
#define EDGEFITTING_H

#include <edgefunction.h>

class edgefitting
{
public:
    edgefitting(int n, BraggEdge::eEdgeFunction ef);
    void intialize_params(double *pars);
    double *get_params();
    void fit(double *x, double *y, int N);
    ~edgefitting();


private:
//    BraggEdge::EdgeFunction edge_lineshape(); /// lineshape to be fitted // does not have to be a class parameter
//    double *x; /// x  axis of the line to be fitted // don't need them either probably
//    double *y; /// y  axis of the line to be fitted
    double *m_pars; /// parameter arrays
    int m_Npars; /// number of parameters
    BraggEdge::eEdgeFunction myfun; /// type of function to be fitted
};

#endif // EDGEFITTING_H
