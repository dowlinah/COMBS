#!/bin/bash

./clean.sh
./build.sh

cd build/bin

(time ./xhpcg) |& tee ./time.txt
valgrind --tool=callgrind ./xhpcg

cd ../..

mv build/bin/callgrind* build/bin/time.txt .
