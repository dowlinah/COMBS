#!/bin/bash
(time ./final_radix_sort_kernel.x 256 10000) |& tee ./time.txt

valgrind --tool=callgrind ./final_radix_sort_kernel.x 256 10000

valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./final_radix_sort_kernel.x 256 10000
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 |& tee mem.txt
