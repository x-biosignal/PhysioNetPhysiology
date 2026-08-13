# Instantaneous phase of a signal

Analytic-signal (Hilbert) instantaneous phase, optionally after
isolating a frequency band. Physiological oscillations (respiration, a
cardiac rhythm) should be reasonably narrow-band for the phase to be
meaningful; pass `band` to bandpass first.

## Usage

``` r
instantaneousPhase(x, sampling_rate = 1, band = NULL)
```

## Arguments

- x:

  Numeric signal.

- sampling_rate:

  Sampling rate in Hz (default 1).

- band:

  Optional `c(low, high)` passband in Hz.

## Value

Numeric vector of wrapped phases in (-pi, pi\].

## Examples

``` r
t <- seq(0, 10, by = 0.01)
ph <- instantaneousPhase(sin(2 * pi * 1 * t), sampling_rate = 100)
```
