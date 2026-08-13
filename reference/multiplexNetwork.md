# Construct a multiplex network

Construct a multiplex network

## Usage

``` r
multiplexNetwork(layers, node_names = NULL)
```

## Arguments

- layers:

  A named list of adjacency matrices, all the same size and node order
  (one layer per coupling channel or frequency band).

- node_names:

  Optional node names (default: the row names of the first layer, or
  `V1..Vn`).

## Value

An object of class `"multiplex"`: a list with `layers`, `nodes`,
`n_layers`, `layer_names`.

## Examples

``` r
L1 <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3)
L2 <- matrix(c(0,0,1, 0,0,0, 1,0,0), 3)
mp <- multiplexNetwork(list(alpha = L1, beta = L2))
```
