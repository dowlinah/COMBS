#!/bin/bash
(time ./sombrero.sh -s small) |& tee ./time.txt

valgrind --tool=callgrind ./sombrero.sh -s small

valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out ./sombrero.sh -s small
grep mem_heap_B massif.out | sed -e 's/mem_heap_B=\(.*\)/\1/' | sort -g | tail -n 1 |& tee mem.txt

mv time.txt callgrind.out* mem.txt ..
