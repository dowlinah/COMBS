#!/bin/bash

cd benchmarks

cd matrix-mpi
make
cd ..
cd phase-retrieval-benchmarks
make
cd ..
cd lid-driven-cavity
make
cd ..
cd FourierBenchmarks
make
cd ..
cd hpcg
./build.sh
cd ../..
