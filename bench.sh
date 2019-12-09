#!/bin/bash

( time $1 ) |& tee ./time.txt
valgrind --tool=callgrind $1
