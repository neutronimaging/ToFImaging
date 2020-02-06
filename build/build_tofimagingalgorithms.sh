#!/bin/bash
if [ `uname` == 'Linux' ]; then
    SPECSTR="-spec linux-g++"
else
    SPECSTR="-spec macx-clang CONFIG+=x86_64"
fi


#dependency: lmfit
DEST_DEPENDENCY=$WORKSPACE/imagingsuite/external/src/lmfit/install
REPOSPATH=$WORKSPACE/ToFImaging
DEST=$WORKSPACE/builds


mkdir $WORKSPACE/imagingsuite/external/src/lmfit/build
mkdir $DEST_DEPENDENCY
cd $WORKSPACE/imagingsuite/external/src/lmfit/build
/usr/bin/cmake -DCMAKE_INSTALL_PREFIX=$DEST_DEPENDENCY .. 
make
make install 
cp $DEST_DEPENDENCY/lib/* $WORKSPACE/lib

mkdir -p $DEST/build-tofimagingalgorithm
cd $DEST/build-tofimagingalgorithm

$QTBINPATH/qmake -makefile -r $SPECSTR -o Makefile ../../ToFImaging/TOF_ImagingAlgorithm/TOF_ImagingAlgorithm.pro

make -f Makefile clean
make -f Makefile mocables all
make -f Makefile

if [ -e "$REPOSPATH/UnitTests" ]
then
    echo "Build tests"

    for f in `ls $REPOSPATH/UnitTests`
    do
        echo "$REPOSPATH/UnitTests/$f/$f.pro"
        if [ -e "$REPOSPATH/UnitTests/$f/$f.pro" ]
        then
            mkdir -p $DEST/build-$f
            cd $DEST/build-$f

            $QTBINPATH/qmake -makefile $SPECSTR -o Makefile ../../ToFImaging/UnitTests/$f/$f.pro
            make -f Makefile clean
            make -f Makefile mocables all
            make -f Makefile
        fi

    done

    echo "Tests built"
fi
