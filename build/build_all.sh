#!/bin/bash

if [ -d "$WORKSPACE/imagingsuite" ] 
then
    echo "Directory $WORKSPACE/imagingsuite exists."
    cd $WORKSPACE/imagingsuite/build/
    ./build_core_kipl.sh
    ./build_core_algorithms.sh
    ./build_lmfit.sh
    cd $WORKSPACE/ToFImaging/build/
else
    echo "Error: Directory $WORKSPACE/ToFImaging does not exists."
fi

./build_tofimagingalgorithms.sh
