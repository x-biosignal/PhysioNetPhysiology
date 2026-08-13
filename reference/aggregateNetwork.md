# Aggregate a multiplex into a single weighted network

Aggregate a multiplex into a single weighted network

## Usage

``` r
aggregateNetwork(x, method = c("sum", "mean"))
```

## Arguments

- x:

  A `"multiplex"` object.

- method:

  `"sum"` (default) or `"mean"` over layers.

## Value

A single weighted adjacency matrix.
