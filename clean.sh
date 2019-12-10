#!/bin/bash

cd benchmarks

cd matrix-mpi
make clean
cd ..
cd phase-retrieval-benchmarks
make clean
cd ..
cd lid-driven-cavity
make clean
cd ..
cd FourierBenchmarks
make clean
cd ..
cd hpcg
./clean.sh
cd 2d-heat
make clean
cd ..
cd sombrero/sombrero-master
make clean
cd ../..
cd fidibench
./build/make clean
./clean.sh
cd ..
cd radix_sort
make clean
cd ../..


