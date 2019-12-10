#!/bin/bash
(time ./build/ctest) |& tee ./time.txt

valgrind --tool=callgrind ./ctest

valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./ctest
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 |& tee mem.txt
