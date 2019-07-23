# Fourier Benchmarks
Benchmark tests the Fastest Fourier Transform in the West (FFTW) Benchmark[1]. This benchmark is critical in the field as solving Fourier Transforms are part of the COMSOL multiphysics simulation software[2]. Since FFTW is an open-source alternative, it would be helpful to compare performance between this library and what is provided in licensed software.

## Installation
Prior to compiling and running the algorithm, some installation steps must be taken. The 'FFTW3' library must be installed as the algorithm computes discrete Fourier transforms. To install the serial library on a Linux machine, enter the following command into a terminal:

sudo apt-get install libfftw3-dev libbftw3-doc

In the FFTW subdirectory, there is a Makefile. Enter ```make``` within the directory and run the output executable.

## Citations
[1] https://github.com/sidneycadot/FourierBenchmarks  
[2] https://www.comsol.com/blogs/how-to-implement-the-fourier-transformation-in-comsol-multiphysics/
