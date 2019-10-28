#include <QtTest>

// add necessary includes here

#include <math/nonlinfit.h>
#include <averageimage.h>
#include <edgefunction.h>
#include <PolynomialCorrection.h>

//#include <sstream>
//#include <iostream>
#include <fstream>


class ToFImagingAlgorithm : public QObject
{
    Q_OBJECT

public:
    ToFImagingAlgorithm();
    ~ToFImagingAlgorithm();

private slots:
    void test_case1();

};

ToFImagingAlgorithm::ToFImagingAlgorithm()
{

}

ToFImagingAlgorithm::~ToFImagingAlgorithm()
{

}

void ToFImagingAlgorithm::test_case1()
{
    short loop=0; //short for loop for input
    string line; //this will contain the data read from the file
    ifstream myfile("/home/carminati_c/git/imagingsuite/core/kipl/UnitTests/data/x.txt"); //opening the file.

    int N=1107;
    double *x = new double[N];
    double *first_guess = new double[N];
    double *computed_first_guess = new double[N];
    double *y = new double[N];

    for (double a; myfile>>a;)
    {
        x[loop]=a;
        loop++;
    }


    ifstream myfile_y ("/home/carminati_c/git/imagingsuite/core/kipl/UnitTests/data/initialmodel.txt"); //opening the file. //path should be related to the lib
    ifstream myfile_y2 ("/home/carminati_c/git/imagingsuite/core/kipl/UnitTests/data/y.txt"); //opening the file. //path should be related to the lib

    short loop_y=0;
    for (double a; myfile_y>>a;)
    {
        first_guess[loop_y]=a;
        loop_y++;

    }

    loop_y=0;
    for (double a; myfile_y2>>a;)
    {
        y[loop_y]=a;
        loop_y++;
    }

    QCOMPARE(loop_y, 1107);

    double *param = new double[7]; // initial parameters
    param[0]=0.056568;
    param[1]=0.0001;
    param[2]=0.0015;
    param[3]=0.315462;
    param[4]=5.3447841;
    param[5]=-0.4700811;
    param[6]=26.929825;

    double *computed_firstedge = new double[N];

//    Nonlinear::Voight vg;

//    ImagingAlgorithms::PolynomialCorrection pc;

    BraggEdge::EdgeFunction myedge(7);

    for (int i=0; i<N; ++i)
    {
        computed_firstedge[i] = BraggEdge::EdgeFunction::EdgeFunctionTExponential(x[i], param);

        }

    // go on here: compare the computed first edge with the loaded one
    // then fitting, and compare the results with the expected final parameters



}

QTEST_APPLESS_MAIN(ToFImagingAlgorithm)

#include "tst_tofimagingalgorithm.moc"
