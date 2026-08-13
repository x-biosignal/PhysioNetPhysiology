# Plot a cardiorespiratory synchrogram

Plot a cardiorespiratory synchrogram

## Usage

``` r
plotSynchrogram(
  x,
  sampling_rate = NULL,
  main = "Cardiorespiratory synchrogram",
  ...
)
```

## Arguments

- x:

  A `"synchrogram"` data frame or a `"cardioresp"` object.

- sampling_rate:

  Sampling rate, to label the x-axis in seconds (optional).

- main:

  Plot title.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`invisible(NULL)`.
