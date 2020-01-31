#ifndef EDGEFITTING_H
#define EDGEFITTING_H

#include <edgefunction.h>
#include <vector>

namespace  ToFImagingAlgorithms

{
class edgefitting
{
public:
    /// \brief class constructor, initialize edge function with n parameters and ef lineshape
    edgefitting(int n, ToFImagingAlgorithms::eEdgeFunction ef);
    /// \brief Initializes fitting parameters, given as input
    void intialize_params(std::vector<double> &pars);
    /// \brief computes initial parameters given the estimated Bragg edge position
    void compute_initial_params(std::vector<double> &x, std::vector<double> &y, int N, double est_t0);
    /// \brief computes intial parameters without the estimated Bragg edge position
    void compute_initial_params(std::vector<double> &x, std::vector<double> &y, int N);
    /// \brief return the current parameters
    void get_params(std::vector<double> &pars);
    /// \brief calls the fitting routine
    void fit(std::vector<double> &x, std::vector<double> &y, int N);
    /// \brief class destructor
    ~edgefitting();


private:
    std::vector<double> m_pars; /// parameter arrays
    int m_Npars; /// number of parameters
    ToFImagingAlgorithms::eEdgeFunction myfun; /// type of function to be fitted
    bool blinear;
};


}


#endif // EDGEFITTING_H
