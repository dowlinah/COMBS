#!/bin/bash

( time ./FFTW/TestFFTW-performance ) |& tee ./time.txt
valgrind --tool=callgrind ./FFTW/TestFFTW-performance
valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./FFTW/TestFFTW-performance
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 | tee ./mem.txt

