#!/bin/bash

( time ./RRR data/data100E 800 .95 .5 100000 100 results100E ) |& tee ./time.txt
valgrind --tool=callgrind ./RRR data/data100E 800 .95 .5 100000 100 results100E 
valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./RRR data/data100E 800 .95 .5 100000 100 results100E 
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 | tee ./mem.txt


