#!/bin/bash

( time mpirun -n 8 ./matrix-mpi 2000 ) |& tee ./time.txt
mpirun -n 8 valgrind --tool=callgrind ./matrix-mpi 2000
