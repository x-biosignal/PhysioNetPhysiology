# State-resolved organ-interaction networks

Builds one TDS network per physiological state and summarises how the
network reconfigures between states. TDS is computed continuously over
the whole recording; each window is assigned to a state, and a pair's
per-state TDS is the percentage of that state's windows that fall inside
a stable coupling period.

## Usage

``` r
tdsNetworkByState(
  X,
  states,
  sampling_rate = 1,
  window = NULL,
  step = NULL,
  max_lag = NULL,
  window_sec = NULL,
  step_sec = NULL,
  max_lag_sec = NULL,
  stability_tol = 1L,
  min_stable = 4L,
  threshold = NULL,
  n_surrogates = 50,
  min_windows = 3L
)
```

## Arguments

- X:

  Node matrix (or any input accepted by
  [`timeDelayStability()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/timeDelayStability.md)).

- states:

  A per-sample state label vector of length `nrow(X)` (character or
  factor). Windows spanning a state boundary take the majority label.

- sampling_rate, window, step, max_lag, window_sec, step_sec,
  max_lag_sec, stability_tol, min_stable:

  Windowing / stability arguments passed to
  [`timeDelayStability()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/timeDelayStability.md).

- threshold:

  Link threshold (TDS percentage). If `NULL` (default), a single global
  surrogate threshold is computed with
  [`tdsSurrogateThreshold()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/tdsSurrogateThreshold.md)
  and applied to every state for comparability.

- n_surrogates:

  Surrogates for the automatic threshold (default 50).

- min_windows:

  Minimum number of windows a state must contain to be analysed (default
  3); smaller states are dropped with a warning.

## Value

An object of class `"tds_state"`: a list with `states` (analysed
levels), `tds` (named list of per-state TDS matrices), `networks` (named
list of
[`tdsNetwork()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/tdsNetwork.md)
results), `threshold`, `window_state`, and `reconfiguration` (a data
frame of per-state link count, density, and mean degree, plus a state x
node degree matrix in `degree`).

## Examples

``` r
set.seed(1)
n <- 4000
x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
half <- n / 2
# heart couples to brain only in state "B"
y <- rnorm(n, sd = 0.3)
y[(half + 6):n] <- x[half:(n - 6)] + rnorm(half - 5, sd = 0.3)
X <- cbind(brain = x, heart = y)
st <- rep(c("A", "B"), each = half)
res <- tdsNetworkByState(X, st, window = 100, step = 50, max_lag = 20,
                         min_stable = 3, threshold = 30)
res$reconfiguration
#>   state windows n_links density mean_degree
#> 1     A      40       0       0           0
#> 2     B      39       1       1           1
```
