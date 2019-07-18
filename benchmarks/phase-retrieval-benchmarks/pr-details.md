# Phase Retrieval Benchmarks
## Description
Benchmarks from the paper, "Benchmark problems for phase retrieval", V. Elser, T.-Y. Lan & T. Bendory.
Used to solve phase retrieval, which uses techniques such as 3D Fourier Transforms
that are applicable to multiphysics simulation designs.
## Installation
Needs the fftw3 library (as does the ks-pde benchmark). Install with the following
line in a terminal window:
```
sudo apt-get install libfftw3-dev libbftw3-doc
```
By default, Ubuntu stores the library in usr/lib. We will need to link the library when compiling. To do so, enter the following in terminal:
```
gcc -O2 RRR.c -lm -L/usr/lib -lfftw3 -o RRR
```
Run with the information given in the benchmark's repository.

## Citations
[1] https://github.com/veitelser/phase-retrieval-benchmarks  
[2] https://arxiv.org/pdf/1706.00399.pdf
[3] https://en.wikipedia.org/wiki/Phase_problem  
[4] http://micro.stanford.edu/wiki/Install_FFTW3  
