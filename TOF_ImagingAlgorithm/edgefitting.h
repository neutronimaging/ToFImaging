#ifndef EDGEFITTING_H
#define EDGEFITTING_H

#include <edgefunction.h>
#include <vector>
#include "tof_imagingalgorithm_global.h"

namespace  ToFImagingAlgorithms

{
class TOF_IMAGINGALGORITHMSHARED_EXPORT edgefitting
{
public:
    /// \brief class constructor, initialize edge function with n parameters and ef lineshape
    edgefitting(int n, ToFImagingAlgorithms::eEdgeFunction ef);
    /// \brief computes initial parameters given the estimated Bragg edge position
    void compute_initial_params(std::vector<double> &x, std::vector<double> &y, int N, double est_t0);
    /// \brief computes intial parameters without the estimated Bragg edge position
    void compute_initial_params(std::vector<double> &x, std::vector<double> &y, int N);
    /// \brief return the current parameters
    void get_params(std::vector<double> &pars);
    /// \brief calls the fitting routine
    void fit(std::vector<double> &x, std::vector<double> &y, int N);
    /// \brief Initializes fitting parameters, given as input
    void initialize_params(std::vector<double> &pars);
    /// \brief calls the fitting routine on selected parameters
    void iterative_fit(std::vector<double> &x, std::vector<double> &y, int N, std::vector<double> &pars);
    /// \brief applies sav-gol filter to smooth edge signal
    void smooth(std::vector<double> &x, std::vector<double> &y);
    /// \brief class destructor
    ~edgefitting();


private:
    std::vector<double> m_pars; ///\param m_pars parameter arrays
    int m_Npars; ///\param m_Npars number of parameters
    ToFImagingAlgorithms::eEdgeFunction myfun; ///\param myfum type of function to be fitted
    bool blinear; ///\param blinear boolean value triggering the computation of linear functions before and after the edge
    bool bsmooth; ///\param bsmooth boolean value triggering the smoothing of the signal, currently not used
    int m_sgWin; ///\param m_sgWin window_lenght for the S-G filter
    int m_sgPoly; ///\param m_sgPoly polynomial order for the S-G filter
};


}


#endif // EDGEFITTING_H
