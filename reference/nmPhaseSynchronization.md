# Scan n:m phase-locking over a grid of ratios

Computes \\\lambda\_{n,m}\\ for every coprime `(n, m)` in the requested
range and returns the dominant ratio. The frequency ratio `f1 : f2`
implied by the winning `(n, m)` is `m : n`.

## Usage

``` r
nmPhaseSynchronization(phi1, phi2, orders_n = 1:8, orders_m = 1:8)
```

## Arguments

- phi1, phi2:

  Instantaneous phases.

- orders_n, orders_m:

  Integer vectors of orders to scan.

## Value

A list with `table` (data frame of `n`, `m`, `lambda`), and the dominant
`n`, `m`, `lambda`, and `ratio` (a `"f1:f2"` string = `m:n`).

## Examples

``` r
t <- seq(0, 120, by = 0.01)
resp <- sin(2 * pi * 0.25 * t)
card <- sin(2 * pi * 1.0 * t)                 # 4x faster: locked 4:1
pr <- instantaneousPhase(resp, 100)
pc <- instantaneousPhase(card, 100)
nmPhaseSynchronization(pc, pr)$ratio
#> [1] "4:1"
```
