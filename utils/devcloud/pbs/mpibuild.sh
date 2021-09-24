#!/bin/bash

if [ $(basename $(pwd)) = "oneSolver" ]
then
    if [[ -d "build" ]]
    then
        echo "Directory build exists. Will be overwritten!"
    else
        echo "Creating build directory ..."
        mkdir "build"
    fi

    cd build
    cmake ..
    cmake --build .
else
   echo "Check current directory. It should be oneSolver or pbs."
fi
