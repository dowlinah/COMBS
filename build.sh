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
cd ..
cd 2d-heat
make
cd ..
cd sombrero/sombrero-master
make
cd ../..
cd fidibench
./build.sh
cd ..
cd radix_sort
make
cd ..
cd monte-carlo
./build.sh
cd ..
cd advection-diffusion
make
cd ..
cd ks-pde
make
cd ../..
