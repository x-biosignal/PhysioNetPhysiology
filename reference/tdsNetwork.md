# Build an organ-interaction network from TDS

Thresholds a TDS matrix into a weighted, undirected network: nodes are
physiological systems and edges are stable dynamic couplings whose
weight is the TDS percentage. Edges below `threshold` are removed.

## Usage

``` r
tdsNetwork(
  x,
  threshold = NULL,
  surrogate = FALSE,
  X = NULL,
  n_surrogates = 50,
  quantile = 0.95
)
```

## Arguments

- x:

  A `"tds"` object from
  [`timeDelayStability()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/timeDelayStability.md),
  or a symmetric numeric TDS matrix.

- threshold:

  Minimum TDS percentage for an edge to be kept. If `NULL` (default) and
  `surrogate = TRUE`, a surrogate threshold is computed with
  [`tdsSurrogateThreshold()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/tdsSurrogateThreshold.md);
  otherwise a value must be supplied.

- surrogate:

  Logical; if `TRUE` and `threshold` is `NULL`, derive the threshold
  from phase-randomised surrogates. Requires `x` to be a `"tds"` object
  (so the underlying series and windowing are available) — supply the
  original node matrix via `X`.

- X:

  Original node matrix, required when `surrogate = TRUE`.

- n_surrogates, quantile:

  Passed to
  [`tdsSurrogateThreshold()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/tdsSurrogateThreshold.md).

## Value

An object of class `"tds_network"`: a list with `adjacency` (weighted,
thresholded), `threshold`, `nodes`, `degree` (number of links per node),
`strength` (summed edge weight per node), and `edges` (a data frame of
the surviving links).

## Examples

``` r
set.seed(1)
n <- 3000
x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
y <- c(rep(0, 5), x[seq_len(n - 5)]) + rnorm(n, sd = 0.3)
z <- rnorm(n)
X <- cbind(brain = x, heart = y, muscle = z)
res <- timeDelayStability(X, window = 100, step = 50, max_lag = 20, min_stable = 3)
net <- tdsNetwork(res, threshold = 40)
net$edges
#>    from    to tds lag
#> 1 brain heart 100   5
```
