#include <QtTest>
#include <QTest>

// add necessary includes here

#include <averageimage.h>
#include <edgefunction.h>
#include <edgefitting.h>
#include <fstream>


class ToFImagingAlgorithm : public QObject
{
    Q_OBJECT

public:
    ToFImagingAlgorithm();
    ~ToFImagingAlgorithm();

private slots:
    void test_TransmissionExp();
    void test_TransmissionLin();
    void test_GradientGaussian();
    void test_AttenuationExp();
    void test_AttenuationLin();

};

ToFImagingAlgorithm::ToFImagingAlgorithm()
{

}

ToFImagingAlgorithm::~ToFImagingAlgorithm()
{

}

void ToFImagingAlgorithm::test_TransmissionExp()
{
    short loop=0; //short for loop for input
    string line; //this will contain the data read from the file
    ifstream myfile("../ToFImaging/UnitTests/test_data/x.txt"); //opening the file.

    unsigned int N=1107;
    double *x = new double[N];
    double *first_guess = new double[N];
    double *y = new double[N];

    double eps=0.0001;

    for (double a; myfile>>a;)
    {
        x[loop]=a;
        loop++;
    }


    ifstream myfile_y ("../ToFImaging/UnitTests/test_data/initialmodel.txt"); //opening the file. //path should be related to the lib
    ifstream myfile_y2 ("../ToFImaging/UnitTests/test_data/y.txt"); //opening the file. //path should be related to the lib

    int loop_y=0;
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

    BraggEdge::EdgeFunction myedge(7);


    for (int i=0; i<N; ++i)
    {
        computed_firstedge[i] = BraggEdge::EdgeFunction::EdgeFunctionTExponential(x[i], param);
//        qDebug() << computed_firstedge[i];
//        qDebug() << first_guess[i];
        QVERIFY(fabs(computed_firstedge[i]-first_guess[i])<eps); // compare the computed first edge with the loaded one, passed

        }


    // then fitting, and compare the results with the expected final parameters

    edgefitting myfit(7, BraggEdge::eEdgeFunction::EdgeTransmissionExponential);
    myfit.intialize_params(param);
    myfit.fit(x,y,N);


    double *updated_params = new double[7];
    updated_params = myfit.get_params();

    // compare with expected output

    double *expected_params = new double[7];
    expected_params[0] = 0.05773708;
    expected_params[1] = 6.1353e-05;
    expected_params[2] = 3.6402e-04;
    expected_params[3] = 0.12457211;
    expected_params[4] = 8.01859477;
    expected_params[5] = 0.08703528;
    expected_params[6] = 17.2455669;

    for (int i=0; i<7; ++i)
    {
        qDebug() << "expected: "  << expected_params[i] << ", computed: " << updated_params[i];
        QVERIFY(fabs(expected_params[i]-updated_params[i])<eps);
    }

    delete [] x;
    delete [] y;
    delete [] expected_params;
    delete [] param;
    delete [] updated_params;
    delete [] first_guess;
    delete [] computed_firstedge;


}

void ToFImagingAlgorithm::test_TransmissionLin()
{

    double *param = new double[7]; // initial parameters
    param[0]=0.056568;
    param[1]=0.0001;
    param[2]=0.0015;
    param[3]=0.6950579360230543;
    param[4]=-2.759274924510319;
    param[5]=0.5660430851212845;
    param[6]=-6.692051420263259;


    double *expected_params = new double[7];
    expected_params[0] = 0.05773708;
    expected_params[1] = 6.2354e-05;
    expected_params[2] = 3.5847e-04;
    expected_params[3] = 0.78972162;
    expected_params[4] = -4.0781558;
    expected_params[5] = 0.48510798;
    expected_params[6] = -5.14692195;


}

void ToFImagingAlgorithm::test_GradientGaussian()
{

}

void ToFImagingAlgorithm::test_AttenuationExp()
{

}

void ToFImagingAlgorithm::test_AttenuationLin()
{

}

QTEST_APPLESS_MAIN(ToFImagingAlgorithm)

#include "tst_tofimagingalgorithm.moc"
