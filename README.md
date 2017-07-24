# SparseQuadratureGrids

[![Build Status](https://travis-ci.org/chriselrod/SparseQuadratureGrids.jl.svg?branch=master)](https://travis-ci.org/chriselrod/SparseQuadratureGrids.jl)

[![Coverage Status](https://coveralls.io/repos/chriselrod/SparseQuadratureGrids.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chriselrod/SparseQuadratureGrids.jl?branch=master)

[![codecov.io](http://codecov.io/github/chriselrod/SparseQuadratureGrids.jl/coverage.svg?branch=master)](http://codecov.io/github/chriselrod/SparseQuadratureGrids.jl?branch=master)

## Introduction

This package implements nested sparse grids, offering both GenzKeister and KronrodPatterson. The univariate rules themselves are saved as JLD files, having been provided by:

> Bourquin, R. Exhaustive search for higher-order Kronrod-Patterson Extensions, 2017. https://github.com/Kronrod-Extensions-Library/kes

Given a current lack of support for hider dimensional sparse arrays, the p-dimensional grids are currently stored as dictionaries. They can be output as a p * n matrix of nodes and n-vector of weights.
The function SplitWeights is the primary interface currently, outputting an object where the nodes and weights have been sorted into positive and negative.

This decision was made to support KernelDensity.jl, so that a combined KDE could be created via a weighted combination of a positive and negative weighted KDE -- potentially useful for plotting purposes.

## Example
