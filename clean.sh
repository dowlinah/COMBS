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
cd ..
cd monte-carlo
./clean.sh
cd ..
cd advection-diffusion
make clean.sh
cd ..
cd ks-pde
make clean
cd ..
cd 2d-heat
make clean
cd ..
cd sombrero/sombrero-master
make clean
./clean.sh
cd ../..
cd fidibench/build
make clean
cd ..
./clean.sh
cd ..
cd radix_sort
make clean
cd ../..


