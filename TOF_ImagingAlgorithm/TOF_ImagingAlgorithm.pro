#-------------------------------------------------
#
# Project created by QtCreator 2019-10-23T11:06:39
#
#-------------------------------------------------

QT       -= gui

TARGET = TOF_ImagingAlgorithm
TEMPLATE = lib
CONFIG += c++11

DEFINES += TOF_IMAGINGALGORITHM_LIBRARY

CONFIG(release, debug|release): DESTDIR    = $$PWD/../../lib
else:CONFIG(debug, debug|release): DESTDIR = $$PWD/../../lib/debug


# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
#        tof_imagingalgorithm.cpp \
    edgefunction.cpp \
    edgefitting.cpp

HEADERS += \
#        tof_imagingalgorithm.h \
        tof_imagingalgorithm_global.h \
    edgefunction.h \
    edgefitting.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}


## Do I need all of this?
unix {
    QMAKE_CXXFLAGS += -fPIC -O2
    unix:!macx {
        QMAKE_CXXFLAGS += -fopenmp
        QMAKE_LFLAGS += -lgomp
        LIBS += -lgomp
    }

    unix:macx {
        INCLUDEPATH += /opt/local/include
        QMAKE_LIBDIR += /opt/local/lib
    }

    LIBS += -lm -lz -L/opt/usr/lib -llmfit
}


CONFIG(release, debug|release)           LIBS += -L$$PWD/../../lib -lkipl -lImagingAlgorithms -lFaddeeva
else:CONFIG(debug, debug|release)        LIBS += -L$$PWD/../../lib/debug/ -lkipl -lImagingAlgorithms -lFaddeeva

INCLUDEPATH += $$PWD/../../imagingsuite/core/kipl/kipl/include $$PWD/../../imagingsuite/core/algorithms/ImagingAlgorithms/include
DEPENDPATH += $$PWD/../../imagingsuite/core/kipl/kipl/include $$PWD/../../imagingsuite/core/algorithms/ImagingAlgorithms/include

INCLUDEPATH += $$PWD/../../imagingsuite/external/src/linalg

INCLUDEPATH += $$PWD/../../imagingsuite/external/src/Fadeeva_erf
DEPENDPATH += $$PWD/../../imagingsuite/external/src/Fadeeva_erf
