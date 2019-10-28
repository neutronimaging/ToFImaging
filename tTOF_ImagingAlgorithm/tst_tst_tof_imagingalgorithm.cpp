
#include <sstream>
#include <iostream>
#include <map>
#include <cmath>

#include <QtTest/QtTest>

// add necessary includes here

#include <lmcurve.h>
#include <fstream>

//#include <math/basicprojector.h>
#include <math/edgefunction.h>

#include <edgefunction.h>
#include <averageimage.h>


#include <base/timage.h>
#include <io/io_fits.h>
#include <io/io_tiff.h>

#include <MorphSpotClean.h>
#include <averageimage.h>
#include <piercingpointestimator.h>
#include <pixelinfo.h>
#include <PolynomialCorrection.h>

#include <edgefunction.h>

using namespace std;
class tst_tof_imagingalgorithm : public QObject
{
    Q_OBJECT

public:
    tst_tof_imagingalgorithm();
    ~tst_tof_imagingalgorithm();

private slots:
    void test_createEdgeLineShape();
    void test_EdgeFitting();

};

tst_tof_imagingalgorithm::tst_tof_imagingalgorithm()
{
    qDebug() << "creating";
}

tst_tof_imagingalgorithm::~tst_tof_imagingalgorithm()
{
    qDebug() << "destroying";
}

void tst_tof_imagingalgorithm::test_createEdgeLineShape()
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

    // compute the initial edge with the implemented model

//    BraggEdge::EdgeFunction my_edge(N);

//    for (int i=0; i<N; ++i)
//    {
//        computed_first_guess[i] =  my_edge.EdgeFunctionTExponential(x[i], param);
//    }

    unsigned long size = 30;
    double *t = new double[size];
//    edgefitting::edgefitting myfit(7);
//    edgefunction::edgefunction myedge(7);

    kipl::math::edgefunction myedge(size, t);

    ImagingAlgorithms::AverageImage avg;




}

void tst_tof_imagingalgorithm::test_EdgeFitting()
{

}

QTEST_APPLESS_MAIN(tst_tof_imagingalgorithm)

#include "tst_tst_tof_imagingalgorithm.moc"
