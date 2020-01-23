QT += testlib
QT -= gui

CONFIG += qt console warn_on depend_includepath testcase
CONFIG -= app_bundle

CONFIG(release, debug|release): DESTDIR    = $$PWD/../../../lib
else:CONFIG(debug, debug|release): DESTDIR = $$PWD/../../../lib/debug

TEMPLATE = app

SOURCES +=  tst_tofimagingalgorithm.cpp

CONFIG(release, debug|release)           LIBS += -L$$PWD/../../../lib -lTOF_ImagingAlgorithm -lkipl -lImagingAlgorithms
else:CONFIG(debug, debug|release)        LIBS += -L$$PWD/../../../lib/debug/ -lTOF_ImagingAlgorithm -lkipl -lImagingAlgorithms

#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib/release/ -lkipl
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib/debug/ -lkipl
#else:unix: LIBS += -L$$PWD/../../../lib/ -lkipl

win32 {
    INCLUDEPATH +=  $$PWD/../../../imagingsuite/external/include/lmfit $$PWD/../../../imagingsuite/external/include
    QMAKE_LIBDIR += $$PWD/../../../imagingsuite/external/lib64/lmfit $$PWD/../../../external/lib64

    LIBS +=  -llmfit

}


INCLUDEPATH += $$PWD/../../../imagingsuite/core/kipl/kipl/include
DEPENDPATH += $$PWD/../../../imagingsuite/core/kipl/kipl/include


#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib/release/ -lImagingAlgorithms
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib/debug/ -lImagingAlgorithms
#else:unix: LIBS += -L$$PWD/../../../lib/ -lImagingAlgorithms

INCLUDEPATH += $$PWD/../../../imagingsuite/core/algorithms/ImagingAlgorithms/include
DEPENDPATH += $$PWD/../../../imagingsuite/core/algorithms/ImagingAlgorithms/include


#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib/release/ -lTOF_ImagingAlgorithm
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib/debug/ -lTOF_ImagingAlgorithm
#else:unix: LIBS += -L$$PWD/../../../lib/ -lTOF_ImagingAlgorithm

INCLUDEPATH += $$PWD/../../TOF_ImagingAlgorithm
DEPENDPATH += $$PWD/../../TOF_ImagingAlgorithm

INCLUDEPATH += $$PWD/../../../imagingsuite/external/src/linalg
