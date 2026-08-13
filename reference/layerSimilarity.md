# Pairwise layer similarity

Cosine similarity between the (upper-triangular) edge-weight vectors of
each pair of layers: 1 = identical coupling structure, 0 = disjoint.

## Usage

``` r
layerSimilarity(x)
```

## Arguments

- x:

  A `"multiplex"` object.

## Value

A `n_layers x n_layers` similarity matrix.
