# Cardiorespiratory synchrogram

The classic Schafer-Rosenblum synchrogram: the respiratory phase
observed at each heartbeat, wrapped into `m` respiratory cycles. Under
n:m locking (n heartbeats per m breaths) the points collapse onto `n`
horizontal bands.

## Usage

``` r
synchrogram(resp_phase, rpeaks, m = 1L)
```

## Arguments

- resp_phase:

  Instantaneous respiratory phase (length = n samples), e.g. from
  [`instantaneousPhase()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/instantaneousPhase.md).

- rpeaks:

  Integer R-peak sample indices.

- m:

  Number of respiratory cycles to wrap into (default 1).

## Value

A data frame with `beat` (R-peak sample index) and `psi` (respiratory
phase at that beat, in `[0, m)`), of class `"synchrogram"`.

## Examples

``` r
t <- seq(0, 120, by = 0.01); fs <- 100
resp_phase <- instantaneousPhase(sin(2 * pi * 0.25 * t), fs)
rpeaks <- which(diff(sign(sin(2 * pi * 1.0 * t))) > 0)   # ~1 Hz beats
sg <- synchrogram(resp_phase, rpeaks, m = 1)
```
