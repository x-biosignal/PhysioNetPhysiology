# Plot the optimal-lag series and stable periods for a pair

Draws the classic Time-Delay-Stability trace: the per-window optimal lag
with windows inside stable periods highlighted. Long flat stretches are
the stable couplings TDS rewards.

## Usage

``` r
plotDelayStability(x, i, j, ...)
```

## Arguments

- x:

  A `"tds"` object from
  [`timeDelayStability()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/timeDelayStability.md).

- i, j:

  The two nodes (names or indices) whose coupling to plot.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`invisible(NULL)`; called for the plot.
