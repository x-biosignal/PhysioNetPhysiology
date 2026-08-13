# Node module roles (Guimera-Amaral)

Combines the within-module degree z-score and participation coefficient
into the seven-role classification of Guimera & Amaral (2005).

## Usage

``` r
moduleRoles(adj, communities)
```

## Arguments

- adj:

  A square adjacency matrix.

- communities:

  A vector of module labels, one per node.

## Value

A data frame with `node`, `z` (within-module degree z-score), `P`
(participation coefficient), and `role`.
