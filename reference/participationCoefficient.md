# Participation coefficient

\\P_i = 1 - \sum_s (\kappa\_{is}/k_i)^2\\, where \\\kappa\_{is}\\ is
node `i`'s summed connection weight to module `s` and \\k_i\\ its total
strength. 0 = all links inside one module; approaches 1 as links spread
evenly over many modules.

## Usage

``` r
participationCoefficient(adj, communities)
```

## Arguments

- adj:

  A square (weighted or binary) adjacency matrix.

- communities:

  A vector of module labels, one per node.

## Value

A named numeric vector of participation coefficients in `[0, 1)`.

## Examples

``` r
adj <- matrix(0, 4, 4); adj[1,2] <- adj[2,1] <- 1; adj[3,4] <- adj[4,3] <- 1
participationCoefficient(adj, c(1, 1, 2, 2))
#> 1 2 3 4 
#> 0 0 0 0 
```
