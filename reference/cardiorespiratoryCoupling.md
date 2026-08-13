# Analyse cardiorespiratory coupling

High-level convenience wrapper: extracts respiratory phase (band-passed
Hilbert) and cardiac phase (from R-peaks), scans n:m phase locking, and
returns the dominant heartbeats-per-breath ratio, its locking index, and
the synchrogram.

## Usage

``` r
cardiorespiratoryCoupling(
  resp,
  rpeaks,
  sampling_rate,
  resp_band = c(0.1, 0.5),
  max_ratio = 8L
)
```

## Arguments

- resp:

  Respiration signal.

- rpeaks:

  Integer R-peak sample indices (same time base as `resp`).

- sampling_rate:

  Sampling rate in Hz.

- resp_band:

  Respiratory passband in Hz (default `c(0.1, 0.5)`).

- max_ratio:

  Largest heartbeats-per-breath ratio to consider (default 8).

## Value

An object of class `"cardioresp"`: a list with `ratio` (a `"heart:resp"`
string), `n`, `m`, `lambda` (locking index), `sync_index` (1:1
phase-locking value), `synchrogram`, and the `table` of scanned ratios.

## Examples

``` r
t <- seq(0, 120, by = 0.01); fs <- 100
resp <- sin(2 * pi * 0.25 * t)
card <- sin(2 * pi * 1.0 * t)                       # 4 beats per breath
rpeaks <- which(diff(sign(card)) > 0)
cr <- cardiorespiratoryCoupling(resp, rpeaks, fs)
cr$ratio
#> [1] "4:1"
```
