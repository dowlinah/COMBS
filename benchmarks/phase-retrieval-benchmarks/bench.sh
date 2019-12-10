#!/bin/bash

( time ./RRR data/data100E 800 .95 .5 1000 5 results100E ) |& tee ./time.txt
valgrind --tool=callgrind ./rrr data/data100e 800 .95 .5 1000 5 results100e 
valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./rrr data/data100e 800 .95 .5 1000 5 results100e 
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 | tee ./mem.txt


