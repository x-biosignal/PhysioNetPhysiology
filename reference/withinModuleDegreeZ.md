# Within-module degree z-score

Standardises each node's within-module strength relative to the other
nodes in its module. High values mark within-module hubs.

## Usage

``` r
withinModuleDegreeZ(adj, communities)
```

## Arguments

- adj:

  A square adjacency matrix.

- communities:

  A vector of module labels, one per node.

## Value

A named numeric vector of z-scores.
