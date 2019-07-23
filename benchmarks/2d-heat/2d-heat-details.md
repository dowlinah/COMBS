# 2D Heat Benchmark
This benchmark compares the utility of MPI, OpenMP, and serial implementations of the 2d-heat benchmark[1]. This work is related to high-performance computing applications as shown in related papers that attempt to parallelize the 2D heat equation to improve performance. The tight coupling of a heat distribution algorithm makes it an obvious candidate for a multiphysics benchmark suite.

## Installation
The src directory includes 3 implementations of the 2D heat algorithm. Compile the desired .c file and run.

## Citations
[1] https://github.com/estrabd/2d-heat  
[2] Horak, V., & Gruber, P. (2005). Parallel Numerical Solution of 2-D Heat Equation. Parallel Numerics, 47-56. Retrieved from https://pdfs.semanticscholar.org/3d15/665a611a4a4f8dc24f77c05025c98cbce3d9.pdf.
