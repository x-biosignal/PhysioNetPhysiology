# Supra-adjacency matrix of a multiplex

Builds the `(N*L) x (N*L)` block matrix with the layer adjacencies on
the diagonal blocks and `interlayer` coupling on the identity
off-diagonal blocks (categorical coupling of a node to its own replicas
across layers).

## Usage

``` r
supraAdjacency(x, interlayer = 1)
```

## Arguments

- x:

  A `"multiplex"` object.

- interlayer:

  Interlayer coupling weight (default 1).

## Value

A `(N*L) x (N*L)` numeric matrix.
