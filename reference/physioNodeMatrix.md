# Assemble a physiological node matrix

Builds the time-by-nodes matrix that
[`timeDelayStability()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/timeDelayStability.md)
consumes, where each node is one physiological variable (an organ system
or a derived feature such as an EEG band-power, heart rate, respiration
amplitude, or an EMG envelope). Inputs are a named list of equal-length
numeric vectors, a numeric matrix, or a data frame.

## Usage

``` r
physioNodeMatrix(
  x,
  nodes = NULL,
  standardize = FALSE,
  assay = NULL,
  channels = NULL,
  feature = c("raw", "envelope"),
  target_rate = NULL
)
```

## Arguments

- x:

  A named list of equal-length numeric vectors, a numeric matrix (time x
  nodes), a data frame of numeric columns, or a `PhysioExperiment` /
  `MultiPhysioExperiment` object (channels/modalities become nodes).

- nodes:

  Optional character vector of node names (defaults to the names of
  `x`).

- standardize:

  If `TRUE`, z-normalise each column (default `FALSE`; TDS z-normalises
  per window internally, so this is usually unnecessary).

- assay:

  For a `PhysioExperiment`, which assay to read (default: its default
  assay).

- channels:

  For a `PhysioExperiment`, an optional subset of channels.

- feature:

  Per-channel feature: `"raw"` (default) or `"envelope"` (the
  analytic-signal amplitude — often the right node series for coupling).

- target_rate:

  Optional rate (Hz) to down-sample each channel to, by block-averaging.
  Required for a `MultiPhysioExperiment` whose experiments have
  differing sampling rates.

## Value

A numeric matrix, time (rows) by nodes (columns), with column names.

## Details

In Network Physiology, coupling is assessed between the *dynamics* of
these variables, so nodes are typically instantaneous features sampled
on a common, relatively slow grid (e.g. 1 Hz), not the raw high-rate
waveforms. Assemble such features with the relevant `Physio*` package
and pass them here.

## Examples

``` r
nodes <- list(
  hr   = cumsum(rnorm(500)),
  resp = sin(seq(0, 20 * pi, length.out = 500)),
  emg  = abs(rnorm(500))
)
M <- physioNodeMatrix(nodes)
dim(M)
#> [1] 500   3
```
