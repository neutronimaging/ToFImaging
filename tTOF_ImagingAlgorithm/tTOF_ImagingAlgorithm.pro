
QT += testlib
QT -= gui

TARGET = tst_tst_tof_imagingalgorithm

CONFIG  +=  console
CONFIG  -=  app_bundle
CONFIG  +=  c++11


CONFIG(release, debug|release): DESTDIR    = $$PWD/../../lib
else:CONFIG(debug, debug|release): DESTDIR = $$PWD/../../lib/debug

TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
#DEFINES += QT_DEPRECATED_WARNINGS
#DEFINES += SRCDIR=\\\"$$PWD/\\\"




#unix {
#    target.path = /usr/lib
#    INSTALLS += target
#}


unix {
    INCLUDEPATH += $$PWD/../../imagingsuite/external/src/linalg
    QMAKE_CXXFLAGS += -fPIC -O2

    unix:!macx {
        QMAKE_CXXFLAGS += -fopenmp
        QMAKE_LFLAGS += -lgomp
        LIBS += -lgomp
    }src

    unix:macx {
        INCLUDEPATH += /opt/local/include
        QMAKE_LIBDIR += /opt/local/lib
    }

    LIBS +=  -lm -lz -ltiff -lfftw3 -lfftw3f -lcfitsio
}

win32 {
    contains(QMAKE_HOST.arch, x86_64):{
    QMAKE_LFLAGS += /MACHINE:X64
    }
    INCLUDEPATH += $$PWD/../../imagingsuite/external/src/linalg $$PWD/../../imagingsuite/external/include $$PWD/../../imagingsuite/external/include/cfitsio
    QMAKE_LIBDIR += $$PWD/../../imagingsuite/external/lib64
    QMAKE_CXXFLAGS += /openmp /O2

     LIBS += -llibtiff -lcfitsio -lzlib_a -llibfftw3-3 -llibfftw3f-3 -lIphlpapi
}

SOURCES += \
        tst_tst_tof_imagingalgorithm.cpp

DEFINES += SRCDIR=\\\"$$PWD/\\\"
CONFIG(debug, debug|release): DEFINES += DEBUG

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0




win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../lib/release/
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../lib/debug/
else:unix: LIBS += -L$$PWD/../../lib/

LIBS += -lkipl -lImagingAlgorithms -lTOF_ImagingAlgorithm -lFaddeeva

INCLUDEPATH += $$PWD/../../imagingsuite/core/kipl/kipl/include
DEPENDPATH += $$PWD/../../imagingsuite/core/kipl/kipl/include

INCLUDEPATH += $$PWD/../../imagingsuite/core/algorithms/ImagingAlgorithms/include
DEPENDPATH += $$PWD/../../imagingsuite/core/algorithms/ImagingAlgorithms/include

INCLUDEPATH += $$PWD/../TOF_ImagingAlgorithm
DEPENDPATH += $$PWD/../TOF_ImagingAlgorithm

INCLUDEPATH += $$PWD/../../imagingsuite/external/src/Fadeeva_erf # to be removed when switching to another imagingsuite branch without Fadeeva
DEPENDPATH += $$PWD/../../imagingsuite/external/src/Fadeeva_erf

#INCLUDEPATH += $$PWD/../../imagingsuite/external/src/linalg


message(test)
message(pwd $$PWD)
message(includepath $$INCLUDEPATH)
message(dependpath $$DEPENDPATH)
message($$PWD/../../lib/)


