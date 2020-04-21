# Benchmark Collection

This repository contains small benchmarks. These benchmarks are mainly used to pre-test smaller code snippets before including them in [`compboost`](https://github.com/schalkdaniel/compboost). The collection contains:

- Fitting a linear model with given Cholesky decomposition vs. inverse matrix: `R/inverse_vs_cholesky.R`
- Fitting a linear model to one categorical feature. This is basically calculating the mean within groups: `R/binary_ridge.R`
