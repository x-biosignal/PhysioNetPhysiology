# Time Delay Stability (TDS) core engine ----------------------------------
# Bashan et al. 2012 (Nat Commun 3:702). For each pair of signals, slide a
# window along the recording; within each window z-normalise both signals and
# find the lag tau0 that maximises |cross-correlation|. A window belongs to a
# "stable" period when tau0 stays within `stability_tol` samples of the previous
# window for a run of at least `min_stable` consecutive windows. TDS strength =
# percentage of windows that fall inside such stable periods.

#' @keywords internal
#' @noRd
.znorm <- function(v) {
  s <- stats::sd(v)
  if (!is.finite(s) || s == 0) return(rep(0, length(v)))
  (v - mean(v)) / s
}

# Windowed peak-lag cross-correlation for one signal pair.
# Returns the per-window optimal lag (tau0), the peak cc value, the window
# start indices, and a logical vector marking windows inside stable periods.
#' @keywords internal
#' @noRd
.pair_tds <- function(x, y, win, step, max_lag,
                      stability_tol = 1L, min_stable = 4L) {
  n <- length(x)
  starts <- seq.int(1L, n - win + 1L, by = step)
  nt <- length(starts)

  pl <- tds_peaklag_cpp(as.double(x), as.double(y),
                        as.integer(starts), as.integer(win),
                        as.integer(max_lag))
  tau0 <- pl$tau0
  peak <- pl$peak

  in_stable <- .stable_windows(tau0, stability_tol, min_stable)
  list(
    tau0 = tau0, peak = peak, starts = starts,
    in_stable = in_stable, n_windows = nt,
    tds = if (nt > 0L) 100 * sum(in_stable) / nt else 0
  )
}

# Pure-R reference for the windowed peak-lag cross-correlation, kept for
# verifying the C++ engine produces identical results.
#' @keywords internal
#' @noRd
.pair_peaklag_r <- function(x, y, starts, win, max_lag) {
  nt <- length(starts)
  lags <- seq.int(-max_lag, max_lag)
  tau0 <- integer(nt); peak <- numeric(nt)
  for (t in seq_len(nt)) {
    idx <- starts[t]:(starts[t] + win - 1L)
    xw <- .znorm(x[idx]); yw <- .znorm(y[idx])
    cc <- numeric(length(lags))
    for (k in seq_along(lags)) {
      lag <- lags[k]
      if (lag >= 0L) { a <- xw[seq_len(win - lag)]; b <- yw[(1L + lag):win] }
      else { a <- xw[(1L - lag):win]; b <- yw[seq_len(win + lag)] }
      cc[k] <- sum(a * b) / win
    }
    j <- which.max(abs(cc))
    tau0[t] <- lags[j]; peak[t] <- cc[j]
  }
  list(tau0 = tau0, peak = peak)
}

# Mark windows belonging to a stable plateau. `step_ok[t]` is TRUE when window t
# is within tolerance of window t-1; a run of L consecutive TRUEs covers L+1
# windows, which counts as stable only if L+1 >= min_stable.
#' @keywords internal
#' @noRd
.stable_windows <- function(tau0, tol = 1L, min_stable = 4L) {
  nt <- length(tau0)
  in_stable <- rep(FALSE, nt)
  if (nt < 2L) return(in_stable)
  step_ok <- c(FALSE, abs(diff(tau0)) <= tol)
  r <- rle(step_ok)
  pos <- 1L
  for (i in seq_along(r$lengths)) {
    L <- r$lengths[i]
    if (isTRUE(r$values[i]) && (L + 1L) >= min_stable) {
      in_stable[(pos - 1L):(pos + L - 1L)] <- TRUE
    }
    pos <- pos + L
  }
  in_stable
}

# Resolve window / step / max_lag from either sample or second units.
#' @keywords internal
#' @noRd
.resolve_windowing <- function(n, sampling_rate,
                               window, step, max_lag,
                               window_sec, step_sec, max_lag_sec) {
  as_samp <- function(x_samp, x_sec, default_sec) {
    if (!is.null(x_samp)) return(as.integer(round(x_samp)))
    if (is.null(x_sec)) x_sec <- default_sec
    max(1L, as.integer(round(x_sec * sampling_rate)))
  }
  win <- as_samp(window, window_sec, 60)
  stp <- as_samp(step, step_sec, NA)
  if (is.null(step) && is.null(step_sec)) stp <- max(1L, win %/% 2L)
  mlag <- as_samp(max_lag, max_lag_sec, NA)
  if (is.null(max_lag) && is.null(max_lag_sec)) mlag <- max(1L, win %/% 4L)

  if (win > n)
    stop(sprintf("window (%d samples) is longer than the series (%d samples).",
                 win, n), call. = FALSE)
  if (mlag >= win)
    stop("max_lag must be smaller than the window length.", call. = FALSE)
  list(window = win, step = stp, max_lag = mlag)
}

#' Time Delay Stability across a multivariate physiological series
#'
#' Computes pairwise Time Delay Stability (TDS) for every pair of columns in
#' `X`. Each column is treated as one physiological node (e.g. an EEG band-power,
#' heart rate, respiration, an EMG envelope). The result is a symmetric matrix
#' of TDS strengths (percentage of time each pair is stably coupled).
#'
#' Windowing may be given in samples (`window`, `step`, `max_lag`) or in seconds
#' (`window_sec`, `step_sec`, `max_lag_sec`, using `sampling_rate`). Sample
#' arguments take precedence. Sensible defaults derive `step` and `max_lag` from
#' the window when unspecified.
#'
#' @param X A numeric matrix or data frame, time (rows) by nodes (columns), or an
#'   object accepted by [physioNodeMatrix()] (named list of equal-length vectors).
#' @param sampling_rate Sampling rate in Hz of the rows of `X` (default 1).
#' @param window,step,max_lag Window length, step, and maximum lag in **samples**
#'   (override the `*_sec` arguments when supplied).
#' @param window_sec,step_sec,max_lag_sec The same quantities in **seconds**.
#'   Defaults: `window_sec = 60`, `step = window/2`, `max_lag = window/4`.
#' @param stability_tol Maximum change in optimal lag (in samples) between
#'   consecutive windows for the pair to be considered stable (default 1).
#' @param min_stable Minimum number of consecutive windows forming a stable
#'   period (default 4).
#' @param nodes Optional character vector of node names (defaults to the column
#'   names of `X`).
#'
#' @return An object of class `"tds"`: a list with `tds` (nodes x nodes TDS
#'   percentage matrix, diagonal 0), `stable_lag` (nodes x nodes dominant stable
#'   lag in samples), `pairs` (per-pair detail incl. `tau0` and `in_stable`),
#'   `nodes`, `n_windows`, and `params`.
#'
#' @examples
#' set.seed(1)
#' n <- 3000
#' x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2))
#' x[is.na(x)] <- 0
#' y <- c(rep(0, 5), x[seq_len(n - 5)]) + rnorm(n, sd = 0.3)  # y lags x by 5
#' z <- rnorm(n)                                              # unrelated
#' X <- cbind(x = x, y = y, z = z)
#' res <- timeDelayStability(X, window = 100, step = 50, max_lag = 20,
#'                           min_stable = 3)
#' res$tds
#'
#' @references Bashan et al. (2012) \doi{10.1038/ncomms1705}.
#' @export
timeDelayStability <- function(X, sampling_rate = 1,
                               window = NULL, step = NULL, max_lag = NULL,
                               window_sec = NULL, step_sec = NULL,
                               max_lag_sec = NULL,
                               stability_tol = 1L, min_stable = 4L,
                               nodes = NULL) {
  X <- .as_node_matrix(X)
  if (ncol(X) < 2L)
    stop("X must have at least two columns (nodes).", call. = FALSE)
  if (anyNA(X))
    stop("X contains NA; impute or drop missing samples first.", call. = FALSE)
  if (is.null(nodes)) nodes <- colnames(X)
  if (is.null(nodes)) nodes <- paste0("V", seq_len(ncol(X)))

  wp <- .resolve_windowing(nrow(X), sampling_rate,
                           window, step, max_lag,
                           window_sec, step_sec, max_lag_sec)
  P <- ncol(X)
  tds_mat <- matrix(0, P, P, dimnames = list(nodes, nodes))
  lag_mat <- matrix(NA_real_, P, P, dimnames = list(nodes, nodes))
  pair_detail <- vector("list", 0L)
  n_windows <- NA_integer_

  for (i in seq_len(P - 1L)) {
    for (j in (i + 1L):P) {
      pr <- .pair_tds(X[, i], X[, j], wp$window, wp$step, wp$max_lag,
                      stability_tol, min_stable)
      tds_mat[i, j] <- tds_mat[j, i] <- pr$tds
      stable_lag <- if (any(pr$in_stable))
        stats::median(pr$tau0[pr$in_stable]) else NA_real_
      lag_mat[i, j] <- stable_lag
      lag_mat[j, i] <- -stable_lag
      pair_detail[[paste(nodes[i], nodes[j], sep = "~")]] <- pr
      n_windows <- pr$n_windows
    }
  }

  structure(
    list(tds = tds_mat, stable_lag = lag_mat, pairs = pair_detail,
         nodes = nodes, n_windows = n_windows,
         params = c(wp, list(sampling_rate = sampling_rate,
                             stability_tol = stability_tol,
                             min_stable = min_stable))),
    class = "tds"
  )
}

#' @export
print.tds <- function(x, ...) {
  cat("<Time Delay Stability>\n")
  cat(sprintf("  nodes:      %d (%s)\n", length(x$nodes),
              paste(utils::head(x$nodes, 6), collapse = ", ")))
  cat(sprintf("  windows:    %d (win=%d, step=%d, max_lag=%d samples)\n",
              x$n_windows, x$params$window, x$params$step, x$params$max_lag))
  off <- x$tds[upper.tri(x$tds)]
  cat(sprintf("  TDS%% range: %.1f - %.1f (mean %.1f)\n",
              min(off), max(off), mean(off)))
  invisible(x)
}

#' @export
summary.tds <- function(object, ...) {
  M <- object$tds
  # index by (row, col) so pair labels and values stay aligned for any n
  idx <- which(upper.tri(M), arr.ind = TRUE)
  nodes <- object$nodes
  df <- data.frame(
    from = nodes[idx[, 1]], to = nodes[idx[, 2]],
    tds = M[idx],
    stable_lag = object$stable_lag[idx],
    row.names = NULL)
  df[order(-df$tds), ]
}
