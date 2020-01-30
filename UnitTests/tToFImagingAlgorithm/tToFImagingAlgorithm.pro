QT += testlib
QT -= gui

CONFIG += qt console warn_on depend_includepath testcase
CONFIG -= app_bundle

TARGET = tst_tofimagingalgorithm

TEMPLATE = app

SOURCES +=  tst_tofimagingalgorithm.cpp

CONFIG += c++11

CONFIG(release, debug|release): DESTDIR    = $$PWD/../../../lib
else:CONFIG(debug, debug|release): DESTDIR = $$PWD/../../../lib/debug

#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib/release/ -lkipl
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib/debug/ -lkipl
#else:unix: LIBS += -L$$PWD/../../../lib/ -lkipl

win32 {

    contains(QMAKE_HOST.arch, x86_64):{
    QMAKE_LFLAGS += /MACHINE:X64
    }

    INCLUDEPATH +=  $$PWD/../../../imagingsuite/external/include/lmfit $$PWD/../../../imagingsuite/external/include
    QMAKE_LIBDIR += $$PWD/../../../imagingsuite/external/lib64/lmfit $$PWD/../../../external/lib64

    LIBS +=  -llmfit

}

#unix:!mac {
#exists(/usr/lib/*NeXus*) {

#    message("-lNeXus exists")
#    DEFINES += HAVE_NEXUS
#    LIBS += -lNeXus -lNeXusCPP
#    SOURCES += ../src/io/io_nexus.cpp
#    HEADERS += ../include/io/io_nexus.h
#}
#else {
#message("-lNeXus does not exists $$HEADERS")
#}

#}

#unix:mac {
#exists($$PWD/../../../imagingsuite/external/mac/lib/*NeXus*) {

#    message("-lNeXus exists")
#    DEFINES += HAVE_NEXUS

#    INCLUDEPATH += $$PWD/../../../imagingsuite/external/mac/include $$PWD/../../../imagingsuite/external/mac/include/nexus $$PWD/../../../imagingsuite/external/mac/include/hdf5
#    DEPENDPATH += $$PWD/../../../imagingsuite/external/mac/include $$PWD/../../../imagingsuite/external/mac/include/nexus $$PWD/../../../imagingsuite/external/mac/include/hdf5

#    LIBS += -L$$PWD/../../../imagingsuite/external/mac/lib -lNeXus.1.0.0 -lNeXusCPP.1.0.0

#}
#else {
#message("-lNeXus does not exists $$HEADERS")
#}

#}

#win32 {

#exists($$PWD/../../../imagingsuite/external/lib64/nexus/*NeXus*) {

#    message("-lNeXus exists")
#    DEFINES += HAVE_NEXUS
#    INCLUDEPATH += $$PWD/../../../imagingsuite/external/include/nexus $$PWD/../../../imagingsuite/external/include/hdf5
#    QMAKE_LIBDIR += $$PWD/../../../imagingsuite/external/lib64/nexus $$PWD/../../../imagingsuite/external/lib64/hdf5

#    LIBS +=  -lNeXus -lNeXusCPP

#    SOURCES += ../src/io/io_nexus.cpp
#    HEADERS += ../include/io/io_nexus.h
#}

#}

CONFIG(release, debug|release)           LIBS += -L$$PWD/../../../lib
else:CONFIG(debug, debug|release)        LIBS += -L$$PWD/../../../lib/debug/

INCLUDEPATH += $$PWD/../../TOF_ImagingAlgorithm
DEPENDPATH += $$PWD/../../TOF_ImagingAlgorithm
LIBS += -lTOF_ImagingAlgorithm

INCLUDEPATH += $$PWD/../../../imagingsuite/core/kipl/kipl/include
DEPENDPATH += $$PWD/../../../imagingsuite/core/kipl/kipl/include
LIBS += -lkipl

INCLUDEPATH += $$PWD/../../../imagingsuite/core/algorithms/ImagingAlgorithms/include
DEPENDPATH += $$PWD/../../../imagingsuite/core/algorithms/ImagingAlgorithms/include
LIBS += -lImagingAlgorithms



INCLUDEPATH += $$PWD/../../../imagingsuite/external/src/linalg
