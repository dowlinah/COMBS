# loops

## Questions

1. Is sum factorization by hand really faster than optimizations by the compiler?
1. Does hand optimization obstruct compiler optimization?
1. How does vectorization play in?
1. Plain loops vs. recursion by dimension?
1. Does it make a difference if we use `u[][][]` or the flattened version `u[]` where we compute indices by hand?
1. How do the answers depend on the number of shape functions/quadrature points?

## Compiler options for `g++`

```
g++ -Ofast -funroll-loops --param max-completely-peeled-insns=10000 --param max-completely-peel-times=10000
```
## Show optimizations and missed optimizations
```
-fopt-info-optimized -fopt-info-missed
```
