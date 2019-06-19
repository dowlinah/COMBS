# KS-PDE: Solving the Kuramoto-Sivashinsky equation
## Description
The Kuramoto-Sivashinsky (KS) equation is one of the principal equations in
connecting partial differential equations (PDEs) and dynamical systems.
Since multiphysics simulations rely heavily on the tight coupling between
dynamic algorithms, this benchmark appears to be reasonable to test
baseline performance, at least on a much simpler, 1D problem. The KS
equation is a highly cited problem in the field dating back over 30
years and being a PDE problem, fits directly in to multiphysics,
since such systems depend on the how quickly a PDE system can
be solved.

## Details/Installing the benchmark
The code for this benchmark is taken from the open-source work done by
Dr. John F. Gibson in an attempt to compare run-time performance between
different programming languages, with the goal of legitimizing Julia as
having the potential to be among the fastest languages, if used correctly.
The benchmark can be adapted to test the performance of just the C code. Prior
to compiling and running the algorithm, some installation steps must be taken.
The 'FFTW3' library must be installed as the algorithm computes discrete
Fourier transforms. To install the serial library on a Linux machine, enter the
following command into a terminal:
```
sudo apt-get install libfftw3-dev libbftw3-doc
```
By default, Ubuntu stores the library in usr/lib. We will need to link the library
when compiling. To do so, enter the following in terminal:
```
clang -O3 -o ksbenchmark-c ksbenchmark.c -L/usr/lib -lfftw3 -lm
```

When running the output file, set Nx as a parameter:
```
./ksbenchmark-c N
```

## Citations
[1]https://github.com/johnfgibson/julia-pde-benchmark  
[2]Hyman, J. M., & Nicolaenko, B. (1986). The Kuramoto-Sivashinsky equation: A bridge between PDES and dynamical systems. Physica D: Nonlinear Phenomena, 18(1-3), 113-126. doi:10.1016/0167-2789(86)90166-1  
[3]http://micro.stanford.edu/wiki/Install_FFTW3
