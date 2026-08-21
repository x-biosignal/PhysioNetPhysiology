# Surrogate-based significance for TDS links ------------------------------
# The TDS of two independent-but-autocorrelated signals is not zero: matching
# spectral content produces spurious stable lags. To get a defensible link
# threshold we destroy the *cross-signal timing* while preserving each signal's
# power spectrum (autocorrelation) via phase randomisation, recompute TDS on the
# surrogate ensemble, and take a high quantile of the surrogate TDS as the
# significance floor.

# Phase-randomised surrogate of a single series: keep the Fourier amplitudes,
# replace the phases with random ones (conjugate-symmetric so the result is
# real). Preserves the linear autocorrelation structure.
#' @keywords internal
#' @noRd
.phase_randomize <- function(v) {
  n <- length(v)
  if (n < 4L) return(sample(v))
  fv <- stats::fft(v)
  amp <- Mod(fv)
  ph <- Arg(fv)
  half <- n %/% 2L
  idx <- 2:half                              # positive freqs excluding DC/Nyquist
  rphase <- stats::runif(length(idx), -pi, pi)
  newph <- ph
  newph[idx] <- rphase
  newph[n - idx + 2L] <- -rphase             # mirror for conjugate symmetry
  # DC (index 1) and, for even n, Nyquist (index half+1) keep their phase (real)
  out <- Re(stats::fft(amp * exp(1i * newph), inverse = TRUE)) / n
  out
}

#' Surrogate significance threshold for TDS links
#'
#' Generates `n_surrogates` phase-randomised surrogate datasets (each column
#' randomised independently, breaking cross-column timing while preserving each
#' column's power spectrum), recomputes the TDS matrix for each, and returns a
#' high quantile of the pooled off-diagonal surrogate TDS values. Links in the
#' real network exceeding this threshold are unlikely to arise from matched
#' spectra alone.
#'
#' @param X Node matrix (or any input accepted by [timeDelayStability()]).
#' @param n_surrogates Number of surrogate datasets (default 50).
#' @param quantile Quantile of the surrogate TDS distribution to use as the
#'   threshold (default 0.95).
#' @param ... Windowing / stability arguments passed to [timeDelayStability()].
#'
#' @return A single numeric TDS percentage threshold, with attribute
#'   `"surrogate_tds"` holding the pooled surrogate values.
#'
#' @examples
#' set.seed(1)
#' n <- 2000
#' x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
#' y <- c(rep(0, 5), x[seq_len(n - 5)]) + rnorm(n, sd = 0.3)
#' X <- cbind(x = x, y = y)
#' thr <- tdsSurrogateThreshold(X, n_surrogates = 20,
#'                              window = 100, step = 50, max_lag = 20,
#'                              min_stable = 3)
#' thr
#'
#' @export
tdsSurrogateThreshold <- function(X, n_surrogates = 50, quantile = 0.95, ...) {
  X <- .as_node_matrix(X)
  vals <- numeric(0)
  for (s in seq_len(n_surrogates)) {
    Xs <- apply(X, 2L, .phase_randomize)
    tds <- timeDelayStability(Xs, ...)$tds
    vals <- c(vals, tds[upper.tri(tds)])
  }
  thr <- as.numeric(stats::quantile(vals, probs = quantile, names = FALSE))
  attr(thr, "surrogate_tds") <- vals
  thr
}
