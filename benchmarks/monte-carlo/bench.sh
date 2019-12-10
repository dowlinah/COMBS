#!/bin/bash

./clean.sh
./build.sh

(time $1 192000000 1000) |& tee ./time.txt
valgrind --tool=callgrind --trace-children=yes $1 192000000 1000
valgrind --tool=massif --trace-children=yes --pages-as-heap=yes --massif-out-file=massif.out $1 192000000 1000
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 | tee ./mem.txt
