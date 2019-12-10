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
./clean.sh
cd ..
cd ks-pde
make clean
cd ../..


