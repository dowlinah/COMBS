# Advection-Diffusion equation
## Description
Advection-diffusion equations are used in the field of computational fluid
dynamics (CFD) to measure phenomena like "the time evolution of chemical or
biological species in a flowing medium such as water or air."
(Hundsdorfer and Verwer). Since such an evolution can be modeled with partial
differential equations (PDEs), an algorithm that models a simpler, 1D aspect
of advection-diffusion can be useful in measuring simulation performance.

## Compiling and running the benchmark
Compile:
```
clang -O3 -o heat-c heat.c -lm
```
Run:
```
./heat-c
```

## Citations
[1] https://github.com/antoine-levitt/benchmark_heat (need to find better benchmark test, still think this is a good benchmark problem)  
[2] Hundsdorfer, W. H., & Verwer, J. G. (2011). Numerical solution of time-dependent advection-diffusion-reaction equations. Berlin: Springer.
