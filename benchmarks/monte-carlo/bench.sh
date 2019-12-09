#!/bin/bash

./clean.sh
./build.sh

(time perl single_core.pl) |& tee ./time.txt
valgrind --tool=callgrind perl single_core.pl
valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out perl single_core.pl
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 | tee ./mem.txt
