# Organ-interaction network from a TDS matrix -----------------------------

#' Build an organ-interaction network from TDS
#'
#' Thresholds a TDS matrix into a weighted, undirected network: nodes are
#' physiological systems and edges are stable dynamic couplings whose weight is
#' the TDS percentage. Edges below `threshold` are removed.
#'
#' @param x A `"tds"` object from [timeDelayStability()], or a symmetric numeric
#'   TDS matrix.
#' @param threshold Minimum TDS percentage for an edge to be kept. If `NULL`
#'   (default) and `surrogate = TRUE`, a surrogate threshold is computed with
#'   [tdsSurrogateThreshold()]; otherwise a value must be supplied.
#' @param surrogate Logical; if `TRUE` and `threshold` is `NULL`, derive the
#'   threshold from phase-randomised surrogates. Requires `x` to be a `"tds"`
#'   object (so the underlying series and windowing are available) — supply the
#'   original node matrix via `X`.
#' @param X Original node matrix, required when `surrogate = TRUE`.
#' @param n_surrogates,quantile Passed to [tdsSurrogateThreshold()].
#'
#' @return An object of class `"tds_network"`: a list with `adjacency` (weighted,
#'   thresholded), `threshold`, `nodes`, `degree` (number of links per node),
#'   `strength` (summed edge weight per node), and `edges` (a data frame of the
#'   surviving links).
#'
#' @examples
#' set.seed(1)
#' n <- 3000
#' x <- as.numeric(stats::filter(rnorm(n), rep(1/5, 5), sides = 2)); x[is.na(x)] <- 0
#' y <- c(rep(0, 5), x[seq_len(n - 5)]) + rnorm(n, sd = 0.3)
#' z <- rnorm(n)
#' X <- cbind(brain = x, heart = y, muscle = z)
#' res <- timeDelayStability(X, window = 100, step = 50, max_lag = 20, min_stable = 3)
#' net <- tdsNetwork(res, threshold = 40)
#' net$edges
#'
#' @export
tdsNetwork <- function(x, threshold = NULL, surrogate = FALSE, X = NULL,
                       n_surrogates = 50, quantile = 0.95) {
  if (inherits(x, "tds")) {
    M <- x$tds
    lag <- x$stable_lag
  } else {
    M <- as.matrix(x)
    if (nrow(M) != ncol(M) || !isTRUE(all.equal(M, t(M))))
      stop("a TDS matrix must be square and symmetric.", call. = FALSE)
    lag <- matrix(NA_real_, nrow(M), ncol(M), dimnames = dimnames(M))
  }
  nodes <- rownames(M) %||% paste0("V", seq_len(nrow(M)))

  if (is.null(threshold)) {
    if (!surrogate)
      stop("supply `threshold`, or set surrogate = TRUE with `X` = node matrix.",
           call. = FALSE)
    if (is.null(X))
      stop("surrogate thresholding needs the original node matrix `X`.",
           call. = FALSE)
    p <- if (inherits(x, "tds")) x$params else list()
    threshold <- tdsSurrogateThreshold(
      X, n_surrogates = n_surrogates, quantile = quantile,
      window = p$window, step = p$step, max_lag = p$max_lag,
      stability_tol = p$stability_tol %||% 1L,
      min_stable = p$min_stable %||% 4L)
  }

  adj <- M
  adj[adj < threshold] <- 0
  diag(adj) <- 0

  ut <- which(upper.tri(adj) & adj > 0, arr.ind = TRUE)
  edges <- if (nrow(ut))
    data.frame(from = nodes[ut[, 1]], to = nodes[ut[, 2]],
               tds = adj[ut], lag = lag[ut], row.names = NULL)
  else
    data.frame(from = character(0), to = character(0),
               tds = numeric(0), lag = numeric(0))
  edges <- edges[order(-edges$tds), ]

  structure(
    list(adjacency = adj, threshold = as.numeric(threshold), nodes = nodes,
         degree = rowSums(adj > 0), strength = rowSums(adj),
         edges = edges, n_links = nrow(edges)),
    class = "tds_network"
  )
}

#' @export
print.tds_network <- function(x, ...) {
  cat("<TDS organ-interaction network>\n")
  cat(sprintf("  nodes:     %d\n", length(x$nodes)))
  cat(sprintf("  threshold: %.1f%% TDS\n", x$threshold))
  cat(sprintf("  links:     %d\n", x$n_links))
  if (x$n_links) {
    top <- utils::head(x$edges, 5)
    cat("  strongest:\n")
    for (i in seq_len(nrow(top)))
      cat(sprintf("    %s -- %s : %.1f%%\n",
                  top$from[i], top$to[i], top$tds[i]))
  }
  invisible(x)
}
