#!/bin/bash

( time $1 ) |& tee ./time.txt
valgrind --tool=callgrind $1
valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out $1
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 | tee ./mem.txt


