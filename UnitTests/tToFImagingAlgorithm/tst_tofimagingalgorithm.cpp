#include <QtTest>
#include <QTest>

// add necessary includes here

#include <averageimage.h>
#include <edgefunction.h>
#include <edgefitting.h>
#include <fstream>
#include <math/gradient.h>
#include <tof2lambda.h>


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
    void test_TOF2lambda();
    void test_lambda2TOF();

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


    ifstream myfile_y ("../ToFImaging/UnitTests/test_data/ini_model_Texp.txt"); //opening the file. //path should be related to the lib
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

    ToFImagingAlgorithms::EdgeFunction myedge(7);


    for (int i=0; i<N; ++i)
    {
        computed_firstedge[i] = ToFImagingAlgorithms::EdgeFunction::EdgeFunctionTExponential(x[i], param);
//        qDebug() << computed_firstedge[i];
//        qDebug() << first_guess[i];
        QVERIFY(fabs(computed_firstedge[i]-first_guess[i])<eps); // compare the computed first edge with the loaded one, passed

        }


    // then fitting, and compare the results with the expected final parameters

   ToFImagingAlgorithms::edgefitting myfit(7, ToFImagingAlgorithms::eEdgeFunction::EdgeTransmissionExponential);
    myfit.intialize_params(param);
    myfit.fit(x,y,N);


    double *updated_params = new double[7];
    myfit.get_params(updated_params);

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

    myfile.close();
    myfile_y.close();
    myfile_y2.close();


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

    ifstream myfile_y ("../ToFImaging/UnitTests/test_data/ini_model_Tlin.txt"); //opening the file. //path should be related to the lib
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




    ToFImagingAlgorithms::EdgeFunction myedge(7);
    double *computed_firstedge = new double[N];
    for (int i=0; i<N; ++i)
    {
        computed_firstedge[i] = ToFImagingAlgorithms::EdgeFunction::EdgeFunctionTLinear(x[i], param);
        QVERIFY(fabs(computed_firstedge[i]-first_guess[i])<eps); // compare the computed first edge with the loaded one, passed

        }

    ToFImagingAlgorithms::edgefitting myfit(7, ToFImagingAlgorithms::eEdgeFunction::EdgeTransmissionLinear);
    myfit.intialize_params(param);
    myfit.fit(x,y,N);


    double *updated_params = new double[7];
    myfit.get_params(updated_params);

    for (int i=0; i<7; ++i)
    {
        qDebug() << "expected: "  << expected_params[i] << ", computed: " << updated_params[i];
        QVERIFY(fabs(expected_params[i]-updated_params[i])<eps);
    }

    myfile.close();
    myfile_y.close();
    myfile_y2.close();


    delete [] x;
    delete [] y;
    delete [] expected_params;
    delete [] param;
    delete [] updated_params;
    delete [] first_guess;
    delete [] computed_firstedge;


}

void ToFImagingAlgorithm::test_GradientGaussian()
{
    // Here I assume that the gradient is already smoothed

    double *param = new double[3];
    param[0] = 0.056568;
    param[1] = 0.0001;
    param[2] = 500.0;

// In case of doubt: put just ones, for the gaussian it works
//    param[0] = 1;
//    param[1] = 1;
//    param[2] = 1;

    double *expected_params = new double[3];
    expected_params[0] = 0.05793304;
    expected_params[1] = 4.6231e-08;
    expected_params[2] = 706.346498;

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

    ifstream myfile_y ("../ToFImaging/UnitTests/test_data/ini_gauss.txt"); //opening the file. //path should be related to the lib
    ifstream myfile_y2 ("../ToFImaging/UnitTests/test_data/smoothed_edge.txt"); //opening the file. //path should be related to the lib


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

    ToFImagingAlgorithms::EdgeFunction myedge(7);
    double *computed_firstedge = new double[N];
    for (int i=0; i<N; ++i)
    {
        computed_firstedge[i] = ToFImagingAlgorithms::EdgeFunction::EdgeGradientGaussian(x[i], param);
        QVERIFY(fabs(computed_firstedge[i]-first_guess[i])<eps); // compare the computed first edge with the loaded one, passed

        }


// check the gradient first

    ifstream myfile_y3 ("../ToFImaging/UnitTests/test_data/gradient.txt");
    double *exp_gradient = new double[N];

    loop_y=0;
    for (double a; myfile_y3>>a;)
    {
        exp_gradient[loop_y]=a;
        loop_y++;

    }

    double *gradient = new double[N];
    kipl::math::num_gradient(y,x,N,gradient);

    for (int i=0; i<N; ++i)
    {
        QVERIFY(fabs(exp_gradient[i]-gradient[i])<eps); // compare the computed first edge with the loaded one, passed

        } // OK. it is the same.




    ToFImagingAlgorithms::edgefitting myfit(3, ToFImagingAlgorithms::eEdgeFunction::EdgeGradientGaussian);
    myfit.intialize_params(param);
    myfit.fit(x,y,N);

    double *updated_params = new double[3];
    myfit.get_params(updated_params);

    for (int i=0; i<3; ++i)
    {
        qDebug() << "expected: "  << expected_params[i] << ", computed: " << updated_params[i];
        QVERIFY(fabs(expected_params[i]-updated_params[i])<eps);
    }

    myfile.close();
    myfile_y.close();
    myfile_y2.close();
    myfile_y3.close();


    delete [] x;
    delete [] y;
    delete [] expected_params;
    delete [] param;
    delete [] updated_params;
    delete [] first_guess;
    delete [] computed_firstedge;




}

void ToFImagingAlgorithm::test_AttenuationExp()
{
    double *param = new double[7];
    double *expected_param = new double[7];

    param[0] = 0.056568;
    param[1] = 0.0001;
    param[2] = 0.0015;
    param[3] = 0.5854773;
    param[4] = 8.245339;
    param[5] = 0.6237067;
    param[6] = -20.09876;

    expected_param[0] = 0.05788515;
    expected_param[1] = 2.1647e-04;
    expected_param[2] = 1.5583e-04;
    expected_param[3] = 0.73414570;
    expected_param[4] = 5.22685127;
    expected_param[5] = 0.44088105;
    expected_param[6] = -16.6111563;

    short loop=0; //short for loop for input
    string line; //this will contain the data read from the file

    ifstream myfile("../ToFImaging/UnitTests/test_data/x.txt"); //opening the file.


    unsigned int N=1107;
    double *x = new double[N];
    double *first_guess = new double[N];
    double *y = new double[N];

    double eps=0.001;

    for (double a; myfile>>a;)
    {
        x[loop]=a;
        loop++;
    }

    ifstream myfile_y ("../ToFImaging/UnitTests/test_data/inimodel_Aexp.txt"); //opening the file. //path should be related to the lib
    ifstream myfile_y2 ("../ToFImaging/UnitTests/test_data/logy.txt"); //opening the file. //path should be related to the lib


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

    ToFImagingAlgorithms::EdgeFunction myedge(7);
    double *computed_firstedge = new double[N];
    for (int i=0; i<N; ++i)
    {
        computed_firstedge[i] = ToFImagingAlgorithms::EdgeFunction::EdgeFunctionAExponential(x[i], param);
//        qDebug() << computed_firstedge[i];
//        qDebug() << first_guess[i];
        QVERIFY(fabs(computed_firstedge[i]-first_guess[i])<eps); // compare the computed first edge with the loaded one, passed

        }


    ToFImagingAlgorithms::edgefitting myfit(7, ToFImagingAlgorithms::eEdgeFunction::EdgeAttenuationExponential);
    myfit.intialize_params(param);
    myfit.fit(x,y,N);


    double *updated_params = new double[7];
    myfit.get_params(updated_params);

    for (int i=0; i<7; ++i)
    {
        qDebug() << "expected: "  << expected_param[i] << ", computed: " << updated_params[i];
//        QVERIFY(fabs(expected_param[i]-updated_params[i])<eps);
    }



}

void ToFImagingAlgorithm::test_AttenuationLin()
{

    double *param = new double[7]; // initial parameters
    param[0]=0.056568;
    param[1]=0.0001;
    param[2]=0.0015;
    param[3]=0.3221873247378136;
    param[4]=5.2674654703300545;
    param[5]=-0.12174127292834003;
    param[6]=31.65784345829804;


    double *expected_params = new double[7];
    expected_params[0] = 0.05788519;
    expected_params[1] = 2.1650e-04;
    expected_params[2] = 1.5581e-04;
    expected_params[3] = 0.14975358;
    expected_params[4] = 7.67236371;
    expected_params[5] = 0.15162227;
    expected_params[6] = 26.4337443;

    short loop=0; //short for loop for input
    string line; //this will contain the data read from the file

    ifstream myfile("../ToFImaging/UnitTests/test_data/x.txt"); //opening the file.


    unsigned int N=1107;
    double *x = new double[N];
    double *first_guess = new double[N];
    double *y = new double[N];

    double eps=0.001;

    for (double a; myfile>>a;)
    {
        x[loop]=a;
        loop++;
    }

    ifstream myfile_y ("../ToFImaging/UnitTests/test_data/inimodel_Alinear.txt"); //opening the file. //path should be related to the lib
    ifstream myfile_y2 ("../ToFImaging/UnitTests/test_data/logy.txt"); //opening the file. //path should be related to the lib


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


    ToFImagingAlgorithms::EdgeFunction myedge(7);
    double *computed_firstedge = new double[N];
    for (int i=0; i<N; ++i)
    {
        computed_firstedge[i] = ToFImagingAlgorithms::EdgeFunction::EdgeFunctionALinear(x[i], param);
//        qDebug() << computed_firstedge[i];
//        qDebug() << first_guess[i];
        QVERIFY(fabs(computed_firstedge[i]-first_guess[i])<eps); // compare the computed first edge with the loaded one, passed

        }

    ToFImagingAlgorithms::edgefitting myfit(7, ToFImagingAlgorithms::eEdgeFunction::EdgeAttenuationLinear);
    myfit.intialize_params(param);
    myfit.fit(x,y,N);


    double *updated_params = new double[7];
    myfit.get_params(updated_params);

    for (int i=0; i<7; ++i)
    {
        qDebug() << "expected: "  << expected_params[i] << ", computed: " << updated_params[i];
//        QVERIFY(fabs(expected_params[i]-updated_params[i])<eps);
    }



}

void ToFImagingAlgorithm::test_TOF2lambda()
{
    // open the tof ifle
    ifstream myfile_tof("../ToFImaging/UnitTests/test_data/tof.txt"); //opening the file.
    ifstream myfile_lambda("../ToFImaging/UnitTests/test_data/lambda.txt"); //opening the file.
    unsigned int N = 300;
    double *exp_lambda = new double[N];
    double *comp_lambda = new double[N];
    double *tof = new double[N];

    int loop=0;
    for (double a; myfile_tof>>a;)
    {
        tof[loop]=a;
        loop++;

    }

    qDebug() << N;
    QCOMPARE(loop, N);

    loop =0;

    for (double a; myfile_lambda>>a;)
    {
        exp_lambda[loop]=a;
        loop++;

    }

    qDebug() << N;
    QCOMPARE(loop, N);

    ToFImagingAlgorithms::ToF2Lambda(tof, comp_lambda, N, 0.0,52);

    for (int i=0; i<N; i++)
    {
        QCOMPARE(comp_lambda[i], exp_lambda[i]);
//        qDebug() << comp_lambda[i];
//        qDebug() << exp_lambda[i];

    }


    delete [] tof;
    delete [] exp_lambda;
    delete [] comp_lambda;
}

void ToFImagingAlgorithm::test_lambda2TOF()
{
    ifstream myfile_tof("../ToFImaging/UnitTests/test_data/tof.txt"); //opening the file.
    ifstream myfile_lambda("../ToFImaging/UnitTests/test_data/lambda.txt"); //opening the file.
    unsigned int N = 300;
    double *exp_tof = new double[N];
    double *comp_tof = new double[N];
    double *lambda = new double[N];

    int loop=0;
    for (double a; myfile_tof>>a;)
    {
        exp_tof[loop]=a;
        loop++;

    }

    qDebug() << N;
    QCOMPARE(loop, N);

    loop =0;

    for (double a; myfile_lambda>>a;)
    {
        lambda[loop]=a;
        loop++;

    }

    qDebug() << N;
    QCOMPARE(loop, N);

    ToFImagingAlgorithms::Lambda2ToF(comp_tof, lambda, N, 0.0,52);

    for (int i=0; i<N; i++)
    {
        QCOMPARE(comp_tof[i], exp_tof[i]);
//        qDebug() << comp_tof[i];
//        qDebug() << exp_tof[i];
    }



    delete [] exp_tof;
    delete [] comp_tof;
    delete [] lambda;

}

QTEST_APPLESS_MAIN(ToFImagingAlgorithm)

#include "tst_tofimagingalgorithm.moc"
