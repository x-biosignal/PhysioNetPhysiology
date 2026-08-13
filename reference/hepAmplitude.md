# Amplitude of the HEP in a component window

Amplitude of the HEP in a component window

## Usage

``` r
hepAmplitude(hep, component = c(0.2, 0.4), fun = c("mean", "peak"))
```

## Arguments

- hep:

  A `"hep"` object.

- component:

  `c(from, to)` seconds defining the component window (default
  `c(0.2, 0.4)`, a common HEP interval).

- fun:

  Summary over the window: `"mean"` (default) or `"peak"` (largest
  absolute value).

## Value

A single amplitude value.
