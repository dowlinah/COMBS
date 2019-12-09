#!/bin/bash
# (time mpirun -n 8 mpi-2dheat.x) |& tee ./time.txt
# (time ./openmp-2dheat.x) |& tee -a ./time.txt
(time ./serial-2dheat.x) |& tee -a ./time.txt

# mpirun -n 8 valgrind --tool=callgrind ./mpi-2dheat.x
# valgrind --tool=callgrind ./openmp-2dheat.x
valgrind --tool=callgrind ./serial-2dheat.x

# mpirun -n 8 valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./mpi-2dheat.x
# grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1
# valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./openmp-2dheat.x
# grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1
valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./serial-2dheat.x
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 |& tee mem.txt