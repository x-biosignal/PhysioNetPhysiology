# Surrogate significance threshold for TDS links

Generates `n_surrogates` phase-randomised surrogate datasets (each
column randomised independently, breaking cross-column timing while
preserving each column's power spectrum), recomputes the TDS matrix for
each, and returns a high quantile of the pooled off-diagonal surrogate
TDS values. Links in the real network exceeding this threshold are
unlikely to arise from matched spectra alone.

## Usage

``` r
tdsSurrogateThreshold(X, n_surrogates = 50, quantile = 0.95, ...)
```

## Arguments

- X:

  Node matrix (or any input accepted by
  [`timeDelayStability()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/timeDelayStability.md)).

- n_surrogates:

  Number of surrogate datasets (default 50).

- quantile:

  Quantile of the surrogate TDS distribution to use as the threshold
  (default 0.95).

- ...:

  Windowing / stability arguments passed to
  [`timeDelayStability()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/timeDelayStability.md).

## Value

A single numeric TDS percentage threshold, with attribute
`"surrogate_tds"` holding the pooled surrogate values.

## Examples

``` r
set.seed(1)
n <- 2000
x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
y <- c(rep(0, 5), x[seq_len(n - 5)]) + rnorm(n, sd = 0.3)
X <- cbind(x = x, y = y)
thr <- tdsSurrogateThreshold(X, n_surrogates = 20,
                             window = 100, step = 50, max_lag = 20,
                             min_stable = 3)
thr
#> [1] 16.15385
#> attr(,"surrogate_tds")
#>  [1]  0.000000 30.769231  7.692308 15.384615  7.692308  0.000000 10.256410
#>  [8]  7.692308  0.000000  0.000000  0.000000  7.692308  0.000000  0.000000
#> [15]  7.692308  7.692308  0.000000  0.000000 10.256410  0.000000
```
