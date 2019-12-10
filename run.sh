#!/bin/bash

cd benchmarks

cd 2d-heat/src
make bench
cd ../..
cd FourierBenchmarks
make bench
cd ..
cd advection-diffusion
make bench
cd ..
cd fidibench
./bench.sh
cd ..
cd hpcg
./bench.sh
cd ..
cd ks-pde
make bench
cd ..
cd lid-driven-cavity
make bench
cd ..
cd matrix-mpi
make bench
cd ..
cd monte-carlo
./bench.sh
cd ..
cd phase-retrieval-benchmarks
make bench
cd ..
cd radix_sort
make bench
cd sombrero/sombrero-master
make bench
cd ../..

cd ../..
./getData.py | tee ./results.txt
