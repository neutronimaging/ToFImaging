#include "tof2lambda.h"
#include <math/mathconstants.h>

namespace ToFImagingAlgorithms {

/// implements TOF to lambda conversion, according to the formula:
/// lambda = h/m (tof-t0)/L/1e-10, where
/// h = 6.62607004e-34 #Planck constant [m^2 kg / s]
/// m = 1.674927471e-27 #Neutron mass [kg]
/// inputs:
/// tof     = pointer to tof array [s]
/// lambda  = pointer to lambda array [A]
/// N       = length of arrays
/// t0      = trigger in the TOF measure [s]
/// L       = flight path [m]
void ToF2Lambda(double *tof, double *lambda, int N, double t0, double L)
{

    for (int i=0; i<N; ++i)
    {
        lambda[i] = dplanck/dnmass*(tof[i]-t0)/(L*1e-10);
    }

}

/// implements lambda to TOF conversion, according to the formula:
/// tof = t0 + (lambda*1e-10)*L*m/h
/// h = 6.62607004e-34 #Planck constant [m^2 kg / s]
/// m = 1.674927471e-27 #Neutron mass [kg]
/// inputs:
/// tof     = pointer to tof array [s]
/// lambda  = pointer to lambda array [A]
/// N       = length of arrays
/// t0      = trigger in the TOF measure [s]
/// L       = flight path [m]
void Lambda2ToF(double *tof, double *lambda, int N, double t0, double L)
{

    for (int i=0; i<N; ++i)
    {
        tof[i]=t0+(lambda[i]*1e-10*L*dnmass)/dplanck;
    }

}

}
