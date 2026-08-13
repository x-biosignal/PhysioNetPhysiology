# n:m phase-locking index

The n:m phase-synchronization index \\\lambda\_{n,m} = \|\langle
e^{i(n\phi_1 - m\phi_2)}\rangle\|\\, in `[0, 1]`. For two rhythms whose
frequency ratio is `f1 : f2 = m : n`, the term `n*phi1 - m*phi2` is
(near-)constant and \\\lambda\\ approaches 1.

## Usage

``` r
phaseSynchronizationIndex(phi1, phi2, n = 1L, m = 1L)
```

## Arguments

- phi1, phi2:

  Instantaneous phases (e.g. from
  [`instantaneousPhase()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/instantaneousPhase.md)).

- n, m:

  Integer orders (default 1:1).

## Value

A single locking index in `[0, 1]`.

## Examples

``` r
t <- seq(0, 60, by = 0.01)
p1 <- instantaneousPhase(sin(2 * pi * 1 * t), 100)
p2 <- instantaneousPhase(sin(2 * pi * 1 * t), 100)
phaseSynchronizationIndex(p1, p2)
#> [1] 1
```
