#!/bin/bash

( time ./FFTW/TestFFTW-performance ) |& tee ./time.txt
valgrind --tool=callgrind ./FFTW/TestFFTW-performance
