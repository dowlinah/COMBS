#!/bin/bash

cd benchmarks

cd matrix-mpi
make bench
cd ..
cd phase-retrieval-benchmarks
make bench
cd ..
cd lid-driven-cavity
make bench
cd ..
cd FourierBenchmarks
make bench
cd ..
cd hpcg
./bench.sh
cd ..
cd 2d-heat/src
make bench
cd ../..
cd sombrero/sombrero-master
make bench
cd ../..
cd fidibench
./bench.sh
cd ..
cd radix_sort
make bench
cd ..
cd monte-carlo
./bench.sh
cd ..
cd advection-diffusion
make bench
cd ..
cd ks-pde
make bench
cd ..
cd moose-next/examples
make test
cd ..
cd ../..
./getData.py | tee ./results.txt
