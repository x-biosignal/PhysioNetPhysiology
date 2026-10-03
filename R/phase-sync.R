# Phase synchronization utilities -----------------------------------------
# Instantaneous phase (analytic signal), n:m phase-locking indices, and the
# scan over frequency ratios. Self-contained (FFT-based Hilbert + brick-wall
# bandpass) so the package carries no heavy signal-processing dependency.

#' @keywords internal
#' @noRd
.analytic <- function(x) {
  n <- length(x)
  X <- stats::fft(x)
  h <- numeric(n)
  if (n %% 2L == 0L) {
    h[1L] <- 1
    h[n %/% 2L + 1L] <- 1
    if (n > 2L) h[2:(n %/% 2L)] <- 2
  } else {
    h[1L] <- 1
    h[2:((n + 1L) %/% 2L)] <- 2
  }
  stats::fft(X * h, inverse = TRUE) / n
}

# Brick-wall FFT bandpass; good enough to isolate a narrow physiological band
# before phase extraction. low/high in Hz; fs sampling rate.
#' @keywords internal
#' @noRd
.fft_bandpass <- function(x, low, high, fs) {
  n <- length(x)
  X <- stats::fft(x)
  k <- 0:(n - 1L)
  k[k > n / 2] <- k[k > n / 2] - n
  freq <- abs(k) * fs / n
  X[freq < low | freq > high] <- 0
  Re(stats::fft(X, inverse = TRUE) / n)
}

#' @keywords internal
#' @noRd
.gcd <- function(a, b) if (b == 0) a else .gcd(b, a %% b)

#' Instantaneous phase of a signal
#'
#' Analytic-signal (Hilbert) instantaneous phase, optionally after isolating a
#' frequency band. Physiological oscillations (respiration, a cardiac rhythm)
#' should be reasonably narrow-band for the phase to be meaningful; pass `band`
#' to bandpass first.
#'
#' @param x Numeric signal.
#' @param sampling_rate Sampling rate in Hz (default 1).
#' @param band Optional `c(low, high)` passband in Hz.
#' @return Numeric vector of wrapped phases in (-pi, pi].
#' @examples
#' t <- seq(0, 10, by = 0.01)
#' ph <- instantaneousPhase(sin(2 * pi * 1 * t), sampling_rate = 100)
#' @export
instantaneousPhase <- function(x, sampling_rate = 1, band = NULL) {
  if (!is.null(band)) {
    stopifnot(length(band) == 2L, band[1] < band[2])
    x <- .fft_bandpass(x, band[1], band[2], sampling_rate)
  }
  Arg(.analytic(x))
}

#' n:m phase-locking index
#'
#' The n:m phase-synchronization index
#' \eqn{\lambda_{n,m} = |\langle e^{i(n\phi_1 - m\phi_2)}\rangle|}, in `[0, 1]`.
#' For two rhythms whose frequency ratio is `f1 : f2 = m : n`, the term
#' `n*phi1 - m*phi2` is (near-)constant and \eqn{\lambda} approaches 1.
#'
#' @param phi1,phi2 Instantaneous phases (e.g. from [instantaneousPhase()]).
#' @param n,m Integer orders (default 1:1).
#' @return A single locking index in `[0, 1]`.
#' @examples
#' t <- seq(0, 60, by = 0.01)
#' p1 <- instantaneousPhase(sin(2 * pi * 1 * t), 100)
#' p2 <- instantaneousPhase(sin(2 * pi * 1 * t), 100)
#' phaseSynchronizationIndex(p1, p2)
#' @export
phaseSynchronizationIndex <- function(phi1, phi2, n = 1L, m = 1L) {
  stopifnot(length(phi1) == length(phi2))
  Mod(mean(exp(1i * (n * phi1 - m * phi2))))
}

#' Scan n:m phase-locking over a grid of ratios
#'
#' Computes \eqn{\lambda_{n,m}} for every coprime `(n, m)` in the requested range
#' and returns the dominant ratio. The frequency ratio `f1 : f2` implied by the
#' winning `(n, m)` is `m : n`.
#'
#' @param phi1,phi2 Instantaneous phases.
#' @param orders_n,orders_m Integer vectors of orders to scan.
#' @return A list with `table` (data frame of `n`, `m`, `lambda`), and the
#'   dominant `n`, `m`, `lambda`, and `ratio` (a `"f1:f2"` string = `m:n`).
#' @examples
#' t <- seq(0, 120, by = 0.01)
#' resp <- sin(2 * pi * 0.25 * t)
#' card <- sin(2 * pi * 1.0 * t)                 # 4x faster: locked 4:1
#' pr <- instantaneousPhase(resp, 100)
#' pc <- instantaneousPhase(card, 100)
#' nmPhaseSynchronization(pc, pr)$ratio
#' @export
nmPhaseSynchronization <- function(phi1, phi2,
                                   orders_n = 1:8, orders_m = 1:8) {
  grid <- expand.grid(n = orders_n, m = orders_m)
  grid <- grid[mapply(function(a, b) .gcd(a, b) == 1L, grid$n, grid$m), ]
  grid$lambda <- mapply(function(n, m)
    phaseSynchronizationIndex(phi1, phi2, n, m), grid$n, grid$m)
  grid <- grid[order(-grid$lambda), ]
  best <- grid[1, ]
  list(table = grid, n = best$n, m = best$m, lambda = best$lambda,
       ratio = sprintf("%d:%d", best$m, best$n))   # f1:f2 = m:n
}
