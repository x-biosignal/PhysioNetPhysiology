# State-resolved TDS networks & reconfiguration ---------------------------
# The signature Network-Physiology result (Bashan et al. 2012): the network of
# organ interactions reorganises across physiological states (wake / light /
# deep / REM sleep). We compute the continuous peak-lag series once over the
# whole recording, attribute each window to the state it falls in, and define
# the per-state TDS of a pair as the fraction of that state's windows that lie
# inside a stable period.

#' @keywords internal
#' @noRd
.majority <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_character_)
  tb <- table(v)
  names(tb)[which.max(tb)]
}

#' State-resolved organ-interaction networks
#'
#' Builds one TDS network per physiological state and summarises how the network
#' reconfigures between states. TDS is computed continuously over the whole
#' recording; each window is assigned to a state, and a pair's per-state TDS is
#' the percentage of that state's windows that fall inside a stable coupling
#' period.
#'
#' @param X Node matrix (or any input accepted by [timeDelayStability()]).
#' @param states A per-sample state label vector of length `nrow(X)` (character
#'   or factor). Windows spanning a state boundary take the majority label.
#' @param sampling_rate,window,step,max_lag,window_sec,step_sec,max_lag_sec,stability_tol,min_stable
#'   Windowing / stability arguments passed to [timeDelayStability()].
#' @param threshold Link threshold (TDS percentage). If `NULL` (default), a
#'   single global surrogate threshold is computed with [tdsSurrogateThreshold()]
#'   and applied to every state for comparability.
#' @param n_surrogates Surrogates for the automatic threshold (default 50).
#' @param min_windows Minimum number of windows a state must contain to be
#'   analysed (default 3); smaller states are dropped with a warning.
#'
#' @return An object of class `"tds_state"`: a list with `states` (analysed
#'   levels), `tds` (named list of per-state TDS matrices), `networks` (named
#'   list of [tdsNetwork()] results), `threshold`, `window_state`, and
#'   `reconfiguration` (a data frame of per-state link count, density, and mean
#'   degree, plus a state x node degree matrix in `degree`).
#'
#' @examples
#' set.seed(1)
#' n <- 4000
#' x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
#' half <- n / 2
#' # heart couples to brain only in state "B"
#' y <- rnorm(n, sd = 0.3)
#' y[(half + 6):n] <- x[half:(n - 6)] + rnorm(half - 5, sd = 0.3)
#' X <- cbind(brain = x, heart = y)
#' st <- rep(c("A", "B"), each = half)
#' res <- tdsNetworkByState(X, st, window = 100, step = 50, max_lag = 20,
#'                          min_stable = 3, threshold = 30)
#' res$reconfiguration
#'
#' @export
tdsNetworkByState <- function(X, states, sampling_rate = 1,
                              window = NULL, step = NULL, max_lag = NULL,
                              window_sec = NULL, step_sec = NULL,
                              max_lag_sec = NULL,
                              stability_tol = 1L, min_stable = 4L,
                              threshold = NULL, n_surrogates = 50,
                              min_windows = 3L) {
  X <- .as_node_matrix(X)
  if (length(states) != nrow(X))
    stop("`states` must have one label per row of X.", call. = FALSE)
  states <- as.factor(states)

  tdsobj <- timeDelayStability(
    X, sampling_rate = sampling_rate,
    window = window, step = step, max_lag = max_lag,
    window_sec = window_sec, step_sec = step_sec, max_lag_sec = max_lag_sec,
    stability_tol = stability_tol, min_stable = min_stable)

  nodes <- tdsobj$nodes
  win <- tdsobj$params$window
  starts <- tdsobj$pairs[[1]]$starts
  nt <- length(starts)

  # majority state per window
  wstate <- factor(
    vapply(starts, function(s) .majority(states[s:(s + win - 1L)]),
           character(1)),
    levels = levels(states))

  # global surrogate threshold if not supplied
  if (is.null(threshold)) {
    p <- tdsobj$params
    threshold <- tdsSurrogateThreshold(
      X, n_surrogates = n_surrogates, quantile = 0.95,
      window = p$window, step = p$step, max_lag = p$max_lag,
      stability_tol = p$stability_tol, min_stable = p$min_stable)
  }

  # per-state TDS matrix from the continuous stability masks
  pair_names <- names(tdsobj$pairs)
  combs <- utils::combn(length(nodes), 2)
  level_counts <- table(wstate)
  keep_levels <- names(level_counts)[level_counts >= min_windows]
  dropped <- setdiff(levels(states), keep_levels)
  if (length(dropped))
    warning(sprintf("states with < %d windows dropped: %s",
                    min_windows, paste(dropped, collapse = ", ")),
            call. = FALSE)

  tds_by_state <- list()
  net_by_state <- list()
  for (lv in keep_levels) {
    in_state <- which(wstate == lv)
    M <- matrix(0, length(nodes), length(nodes),
                dimnames = list(nodes, nodes))
    for (c in seq_len(ncol(combs))) {
      i <- combs[1, c]; j <- combs[2, c]
      pr <- tdsobj$pairs[[pair_names[c]]]
      denom <- length(in_state)
      val <- if (denom > 0) 100 * sum(pr$in_stable[in_state]) / denom else 0
      M[i, j] <- M[j, i] <- val
    }
    tds_by_state[[lv]] <- M
    net_by_state[[lv]] <- tdsNetwork(M, threshold = threshold)
  }

  # reconfiguration summary
  nn <- length(nodes)
  max_links <- nn * (nn - 1) / 2
  recon <- data.frame(
    state = keep_levels,
    windows = as.integer(level_counts[keep_levels]),
    n_links = vapply(net_by_state, function(z) z$n_links, integer(1)),
    density = vapply(net_by_state,
                     function(z) z$n_links / max_links, numeric(1)),
    mean_degree = vapply(net_by_state,
                         function(z) mean(z$degree), numeric(1)),
    row.names = NULL)
  degmat <- t(vapply(net_by_state, function(z) z$degree, numeric(nn)))
  rownames(degmat) <- keep_levels

  structure(
    list(states = keep_levels, tds = tds_by_state, networks = net_by_state,
         threshold = as.numeric(threshold), window_state = wstate,
         reconfiguration = recon, degree = degmat, nodes = nodes),
    class = "tds_state")
}

#' @export
print.tds_state <- function(x, ...) {
  cat("<State-resolved TDS networks>\n")
  cat(sprintf("  states:    %s\n", paste(x$states, collapse = ", ")))
  cat(sprintf("  threshold: %.1f%% TDS\n", x$threshold))
  cat("  links per state:\n")
  for (i in seq_len(nrow(x$reconfiguration))) {
    r <- x$reconfiguration[i, ]
    cat(sprintf("    %-8s %d links (density %.2f)\n",
                r$state, r$n_links, r$density))
  }
  invisible(x)
}

#' Link-level reconfiguration between two states
#'
#' Returns, for every node pair, the TDS in each of two states and their
#' difference — the edges that appear, vanish, or strengthen between states.
#'
#' @param x A `"tds_state"` object.
#' @param from,to State labels to compare.
#'
#' @return A data frame of `from`, `to` (node pair), the TDS in each state, and
#'   `delta`, ordered by absolute change.
#' @export
tdsReconfiguration <- function(x, from, to) {
  stopifnot(inherits(x, "tds_state"))
  if (!all(c(from, to) %in% x$states))
    stop("`from`/`to` must be analysed states.", call. = FALSE)
  A <- x$tds[[from]]; B <- x$tds[[to]]
  nodes <- x$nodes
  combs <- utils::combn(length(nodes), 2)
  df <- data.frame(
    from_node = nodes[combs[1, ]], to_node = nodes[combs[2, ]],
    tds_from = A[t(combs)], tds_to = B[t(combs)], row.names = NULL)
  df$delta <- df$tds_to - df$tds_from
  names(df)[3:4] <- c(paste0("tds_", from), paste0("tds_", to))
  df[order(-abs(df$delta)), ]
}
