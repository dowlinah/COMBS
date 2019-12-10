#!/bin/bash

( time mpirun -n 8 ./matrix-mpi 2000 ) |& tee ./time.txt
mpirun -n 8 valgrind --tool=callgrind ./matrix-mpi 2000
mpirun -n 8 valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./matrix-mpi 2000
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 | tee ./mem.txt


