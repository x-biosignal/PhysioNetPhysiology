# Visualisation -----------------------------------------------------------

#' Plot the optimal-lag series and stable periods for a pair
#'
#' Draws the classic Time-Delay-Stability trace: the per-window optimal lag with
#' windows inside stable periods highlighted. Long flat stretches are the stable
#' couplings TDS rewards.
#'
#' @param x A `"tds"` object from [timeDelayStability()].
#' @param i,j The two nodes (names or indices) whose coupling to plot.
#' @param ... Passed to [graphics::plot()].
#' @return `invisible(NULL)`; called for the plot.
#' @export
plotDelayStability <- function(x, i, j, ...) {
  stopifnot(inherits(x, "tds"))
  key <- .pair_key(x, i, j)
  pr <- x$pairs[[key]]
  if (is.null(pr)) stop("no such node pair.", call. = FALSE)
  w <- seq_along(pr$tau0)
  graphics::plot(w, pr$tau0, type = "l", col = "grey60",
                 xlab = "window", ylab = "optimal lag (samples)",
                 main = sprintf("TDS %s (%.1f%% stable)", gsub("~", " - ", key),
                                pr$tds), ...)
  graphics::points(w[pr$in_stable], pr$tau0[pr$in_stable],
                   pch = 16, col = "firebrick", cex = 0.7)
  graphics::legend("topright", bty = "n", pch = c(NA, 16),
                   lty = c(1, NA), col = c("grey60", "firebrick"),
                   legend = c("optimal lag", "stable"))
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.pair_key <- function(x, i, j) {
  nm <- x$nodes
  a <- if (is.character(i)) i else nm[i]
  b <- if (is.character(j)) j else nm[j]
  k1 <- paste(a, b, sep = "~"); k2 <- paste(b, a, sep = "~")
  if (!is.null(x$pairs[[k1]])) k1 else k2
}

#' Plot a TDS organ-interaction network
#'
#' Circular node layout; edge width and colour encode TDS strength.
#'
#' @param x A `"tds_network"` object.
#' @param main Plot title.
#' @param ... Ignored.
#' @return `invisible(NULL)`.
#' @export
plotTDSnetwork <- function(x, main = "TDS organ-interaction network", ...) {
  stopifnot(inherits(x, "tds_network"))
  nodes <- x$nodes; nn <- length(nodes)
  ang <- seq(0, 2 * pi, length.out = nn + 1)[seq_len(nn)]
  px <- cos(ang); py <- sin(ang)
  op <- graphics::par(mar = c(1, 1, 2, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3),
                 asp = 1, axes = FALSE, xlab = "", ylab = "", main = main)
  adj <- x$adjacency
  if (x$n_links) {
    wmax <- max(adj)
    for (a in seq_len(nn - 1)) for (b in (a + 1):nn) {
      w <- adj[a, b]
      if (w > 0) {
        graphics::segments(px[a], py[a], px[b], py[b],
                           lwd = 1 + 4 * w / wmax,
                           col = grDevices::adjustcolor("steelblue",
                                                        0.3 + 0.7 * w / wmax))
      }
    }
  }
  graphics::points(px, py, pch = 21, bg = "white", cex = 3, lwd = 2,
                   col = "grey30")
  graphics::text(px * 1.18, py * 1.18, labels = nodes, cex = 0.9)
  invisible(NULL)
}

#' Plot network reconfiguration across states
#'
#' Heatmap of per-state node degree (states x nodes): how each system's number
#' of stable couplings changes across physiological states.
#'
#' @param x A `"tds_state"` object.
#' @param main Plot title.
#' @param ... Ignored.
#' @return `invisible(NULL)`.
#' @export
plotTDSreconfiguration <- function(x, main = "TDS network reconfiguration", ...) {
  stopifnot(inherits(x, "tds_state"))
  D <- x$degree
  ns <- nrow(D); nn <- ncol(D)
  op <- graphics::par(mar = c(4, 6, 3, 2)); on.exit(graphics::par(op))
  pal <- grDevices::colorRampPalette(c("#f7fbff", "#08519c"))(32)
  graphics::image(x = seq_len(nn), y = seq_len(ns), z = t(D),
                  col = pal, axes = FALSE, xlab = "", ylab = "",
                  main = main, zlim = c(0, max(D)))
  graphics::axis(1, at = seq_len(nn), labels = colnames(D), las = 2, cex.axis = 0.8)
  graphics::axis(2, at = seq_len(ns), labels = rownames(D), las = 1, cex.axis = 0.9)
  graphics::box()
  for (r in seq_len(ns)) for (cc in seq_len(nn))
    graphics::text(cc, r, D[r, cc], cex = 0.8,
                   col = if (D[r, cc] > max(D) / 2) "white" else "grey20")
  invisible(NULL)
}
