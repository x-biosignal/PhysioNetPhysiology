# Time Delay Stability across a multivariate physiological series

Computes pairwise Time Delay Stability (TDS) for every pair of columns
in `X`. Each column is treated as one physiological node (e.g. an EEG
band-power, heart rate, respiration, an EMG envelope). The result is a
symmetric matrix of TDS strengths (percentage of time each pair is
stably coupled).

## Usage

``` r
timeDelayStability(
  X,
  sampling_rate = 1,
  window = NULL,
  step = NULL,
  max_lag = NULL,
  window_sec = NULL,
  step_sec = NULL,
  max_lag_sec = NULL,
  stability_tol = 1L,
  min_stable = 4L,
  nodes = NULL
)
```

## Arguments

- X:

  A numeric matrix or data frame, time (rows) by nodes (columns), or an
  object accepted by
  [`physioNodeMatrix()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/physioNodeMatrix.md)
  (named list of equal-length vectors).

- sampling_rate:

  Sampling rate in Hz of the rows of `X` (default 1).

- window, step, max_lag:

  Window length, step, and maximum lag in **samples** (override the
  `*_sec` arguments when supplied).

- window_sec, step_sec, max_lag_sec:

  The same quantities in **seconds**. Defaults: `window_sec = 60`,
  `step = window/2`, `max_lag = window/4`.

- stability_tol:

  Maximum change in optimal lag (in samples) between consecutive windows
  for the pair to be considered stable (default 1).

- min_stable:

  Minimum number of consecutive windows forming a stable period (default
  4).

- nodes:

  Optional character vector of node names (defaults to the column names
  of `X`).

## Value

An object of class `"tds"`: a list with `tds` (nodes x nodes TDS
percentage matrix, diagonal 0), `stable_lag` (nodes x nodes dominant
stable lag in samples), `pairs` (per-pair detail incl. `tau0` and
`in_stable`), `nodes`, `n_windows`, and `params`.

## Details

Windowing may be given in samples (`window`, `step`, `max_lag`) or in
seconds (`window_sec`, `step_sec`, `max_lag_sec`, using
`sampling_rate`). Sample arguments take precedence. Sensible defaults
derive `step` and `max_lag` from the window when unspecified.

## References

Bashan et al. (2012)
[doi:10.1038/ncomms1705](https://doi.org/10.1038/ncomms1705) .

## Examples

``` r
set.seed(1)
n <- 3000
x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2))
x[is.na(x)] <- 0
y <- c(rep(0, 5), x[seq_len(n - 5)]) + rnorm(n, sd = 0.3)  # y lags x by 5
z <- rnorm(n)                                              # unrelated
X <- cbind(x = x, y = y, z = z)
res <- timeDelayStability(X, window = 100, step = 50, max_lag = 20,
                          min_stable = 3)
res$tds
#>           x          y         z
#> x   0.00000 100.000000 15.254237
#> y 100.00000   0.000000  5.084746
#> z  15.25424   5.084746  0.000000
```
