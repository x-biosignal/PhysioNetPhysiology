# Multiplex participation coefficient

How evenly a node's connections spread across layers. 0 = all links
confined to one layer; 1 = links distributed uniformly across all layers
(Battiston et al. 2014).

## Usage

``` r
multiplexParticipation(x)
```

## Arguments

- x:

  A `"multiplex"` object.

## Value

A named numeric vector per node, in `[0, 1]`.
