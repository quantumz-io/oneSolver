#!/bin/bash

INFILE=benchmarks/exhaustive_search/examples/csp13.qubo
OUTFILE=csp13.csv
MPIPROC=4
WITHGPU=false

if [ "$WITHGPU" = true ]; then
    mpirun -n $MPIPROC ./build/src/one-solver-exhaustive --device-type gpu --input $INFILE --output $OUTFILE
elif [ "$WITHGPU" = false ]; then
    mpirun -n $MPIPROC ./build/src/one-solver-exhaustive --input $INFILE --output $OUTFILE
else
    echo "Unknown device option"
fi
