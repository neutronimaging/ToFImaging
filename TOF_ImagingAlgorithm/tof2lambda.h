#ifndef TOF2LAMBDA_H
#define TOF2LAMBDA_H


#include <type_traits>


namespace ToFImagingAlgorithms {


///\brief performs TOF [s] to lambda [A] conversion
void ToF2Lambda(double *tof, double *lambda, int N, double t0, double L);

///\brief performs lambda [A] to TOF [s] conversion
void Lambda2ToF(double *tof, double *lambda, int N,  double t0, double L);

// Think about the templated versions

/////\brief performs TOF [s] to lambda [A] conversion
//template <typename myType>
//void ToF2Lambda(myType *tof, myType *lambda, int N, myType t0, myType L)
//{
//    static_assert(std::is_floating_point<myType>::value,
//                    "ToF2Lambda function can only be instantiated with floating point types");


//}

/////\brief performs lambda [A] to TOF [s] conversion
//template <typename myType>
//void Lambda2ToF(myType *tof, myType *lambda, int N,  myType t0, myType L)
//{
//    static_assert(std::is_floating_point<myType>::value,
//                    "Lambda2ToF function can only be instantiated with floating point types");

//}

}


#endif // TOF2LAMBDA_H
