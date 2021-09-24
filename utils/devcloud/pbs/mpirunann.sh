#!/bin/bash

infile=benchmarks/exhaustive_search/examples/csp13.qubo
outfile=csp13.csv
MPIPROC=4
WITHGPU=false

if [ "$WITHGPU" = true ]; then
    mpirun -n $MPIPROC ./build/src/one-solver-anneal --device-type gpu --input $infile --output $outfile
elif [ "$WITHGPU" = false ]; then
    mpirun -n $MPIPROC ./build/src/one-solver-anneal --input $infile --output $outfile
else
    echo "Unknown device option"
fi
