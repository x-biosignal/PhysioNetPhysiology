# Permutation significance of the HEP component

Tests whether the HEP component amplitude is larger than expected if the
brain signal were not time-locked to the heartbeat, by recomputing the
component amplitude for `n_perm` datasets with circularly shifted
R-peaks.

## Usage

``` r
hepSignificance(
  eeg,
  rpeaks,
  sampling_rate,
  window = c(-0.1, 0.6),
  baseline = c(-0.1, 0),
  component = c(0.2, 0.4),
  n_perm = 200
)
```

## Arguments

- eeg, rpeaks, sampling_rate, window, baseline:

  As in
  [`heartbeatEvokedPotential()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/heartbeatEvokedPotential.md).

- component:

  Component window in seconds (see
  [`hepAmplitude()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/hepAmplitude.md)).

- n_perm:

  Number of permutations (default 200).

## Value

A list with `observed` (\|component amplitude\|), `null` (permutation
amplitudes), and `p_value` (one-sided).
