#include <QtTest>
//#include <QDebug>

// add necessary includes here

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

}

void tst_tof_imagingalgorithm::test_EdgeFitting()
{

}

QTEST_APPLESS_MAIN(tst_tof_imagingalgorithm)

#include "tst_tst_tof_imagingalgorithm.moc"
