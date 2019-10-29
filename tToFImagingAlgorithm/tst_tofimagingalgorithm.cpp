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

    double eps=0.0001;

    for (double a; myfile>>a;)
    {
        x[loop]=a;
        loop++;
    }


    ifstream myfile_y ("/home/carminati_c/git/imagingsuite/core/kipl/UnitTests/data/initialmodel.txt"); //opening the file. //path should be related to the lib
    ifstream myfile_y2 ("/home/carminati_c/git/imagingsuite/core/kipl/UnitTests/data/y.txt"); //opening the file. //path should be related to the lib

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
        QVERIFY(fabs(computed_firstedge[i]-first_guess[i])<eps); // compare the computed first edge with the loaded one, passed

        }


    // then fitting, and compare the results with the expected final parameters

    edgefitting myfit(7);
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


}

QTEST_APPLESS_MAIN(ToFImagingAlgorithm)

#include "tst_tofimagingalgorithm.moc"
