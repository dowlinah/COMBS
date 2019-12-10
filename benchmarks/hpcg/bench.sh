#!/bin/bash

./clean.sh
./build.sh

cd build/bin

(time ./xhpcg) |& tee ./time.txt
valgrind --tool=callgrind ./xhpcg
valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./xhpcg
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 | tee ./mem.txt

cd ../..

mv build/bin/callgrind* build/bin/time.txt build/bin/mem.txt build/bin/massif.out .
