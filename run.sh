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
cd 2d-heat
make bench
cd ..
cd sombrero/sombrero-master
make bench
cd ../..
cd fidibench
make bench
cd ..
cd radix_sort
make bench
cd ..
cd ../..

./getData.py | tee ./results.txt
