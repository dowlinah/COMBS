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
cd ../..

./getData.py | tee ./results.txt
