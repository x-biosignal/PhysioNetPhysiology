# Link-level reconfiguration between two states

Returns, for every node pair, the TDS in each of two states and their
difference — the edges that appear, vanish, or strengthen between
states.

## Usage

``` r
tdsReconfiguration(x, from, to)
```

## Arguments

- x:

  A `"tds_state"` object.

- from, to:

  State labels to compare.

## Value

A data frame of `from`, `to` (node pair), the TDS in each state, and
`delta`, ordered by absolute change.
