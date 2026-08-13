# Cardiac phase from R-peak marks

Builds a continuous cardiac phase that advances linearly from 0 to
\\2\pi\\ between consecutive R-peaks (wrapped per beat). Samples outside
the first/last beat are `NA`.

## Usage

``` r
cardiacPhaseFromPeaks(rpeaks, n)
```

## Arguments

- rpeaks:

  Integer sample indices of R-peaks (1-based).

- n:

  Length of the output phase vector (number of samples).

## Value

Numeric vector of wrapped cardiac phase in `[0, 2*pi)`, length `n`.

## Examples

``` r
rpeaks <- seq(10, 1000, by = 25)
ph <- cardiacPhaseFromPeaks(rpeaks, 1000)
```
