# Heartbeat-evoked potential (HEP)

Epochs a brain signal around each R-peak, baseline-corrects, and
averages to the heartbeat-evoked potential.

## Usage

``` r
heartbeatEvokedPotential(
  eeg,
  rpeaks,
  sampling_rate,
  window = c(-0.1, 0.6),
  baseline = c(-0.1, 0)
)
```

## Arguments

- eeg:

  Numeric brain signal (one channel).

- rpeaks:

  Integer R-peak sample indices (same time base as `eeg`).

- sampling_rate:

  Sampling rate in Hz.

- window:

  Epoch window `c(pre, post)` in seconds relative to the R-peak (default
  `c(-0.1, 0.6)`).

- baseline:

  Baseline window `c(from, to)` in seconds for per-epoch correction;
  `NULL` to skip (default `c(-0.1, 0)`).

## Value

An object of class `"hep"`: a list with `time` (s, relative to the
R-peak), `hep` (averaged waveform), `epochs` (epochs x time matrix), `n`
(number of epochs), and `sampling_rate`.

## Examples

``` r
set.seed(1); fs <- 250; n <- 30 * fs
rpeaks <- seq(fs, n - fs, by = round(0.9 * fs))
eeg <- rnorm(n)
hep <- heartbeatEvokedPotential(eeg, rpeaks, fs)
```
