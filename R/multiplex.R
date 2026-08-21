# Multiplex / multilayer physiological networks ---------------------------
# Organ interactions are not single-layer: brain waves (delta..beta), heart,
# respiration and muscle couple through several channels at once. A multiplex
# network stacks one adjacency layer per channel over a shared node set
# (Bartsch et al. 2015). This module builds the stack and computes the metrics
# that only make sense across layers.

#' Construct a multiplex network
#'
#' @param layers A named list of adjacency matrices, all the same size and node
#'   order (one layer per coupling channel or frequency band).
#' @param node_names Optional node names (default: the row names of the first
#'   layer, or `V1..Vn`).
#' @return An object of class `"multiplex"`: a list with `layers`, `nodes`,
#'   `n_layers`, `layer_names`.
#' @examples
#' L1 <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3)
#' L2 <- matrix(c(0,0,1, 0,0,0, 1,0,0), 3)
#' mp <- multiplexNetwork(list(alpha = L1, beta = L2))
#' @export
multiplexNetwork <- function(layers, node_names = NULL) {
  if (!is.list(layers) || length(layers) < 2L)
    stop("`layers` must be a list of at least two adjacency matrices.",
         call. = FALSE)
  layers <- lapply(layers, as.matrix)
  d <- dim(layers[[1]])
  if (d[1] != d[2]) stop("layers must be square.", call. = FALSE)
  if (!all(vapply(layers, function(m) all(dim(m) == d), logical(1))))
    stop("all layers must have the same dimensions.", call. = FALSE)
  nodes <- node_names %||% rownames(layers[[1]]) %||% paste0("V", seq_len(d[1]))
  lnames <- names(layers) %||% paste0("layer", seq_along(layers))
  layers <- lapply(layers, function(m) { dimnames(m) <- list(nodes, nodes); m })
  structure(list(layers = layers, nodes = nodes, n_layers = length(layers),
                 layer_names = lnames), class = "multiplex")
}

#' Per-layer and overlapping (multiplex) degree
#'
#' @param x A `"multiplex"` object.
#' @return A list with `per_layer` (nodes x layers degree matrix) and
#'   `overlapping` (summed degree across layers per node).
#' @export
multiplexDegree <- function(x) {
  stopifnot(inherits(x, "multiplex"))
  per <- vapply(x$layers, function(m) rowSums(m > 0), numeric(length(x$nodes)))
  colnames(per) <- x$layer_names
  list(per_layer = per, overlapping = rowSums(per))
}

#' Multiplex participation coefficient
#'
#' How evenly a node's connections spread across layers. 0 = all links confined
#' to one layer; 1 = links distributed uniformly across all layers
#' (Battiston et al. 2014).
#'
#' @param x A `"multiplex"` object.
#' @return A named numeric vector per node, in `[0, 1]`.
#' @export
multiplexParticipation <- function(x) {
  stopifnot(inherits(x, "multiplex"))
  L <- x$n_layers
  per <- multiplexDegree(x)$per_layer          # nodes x layers
  oi <- rowSums(per)
  P <- (L / (L - 1)) * (1 - rowSums((per / oi)^2))
  P[oi == 0] <- 0
  stats::setNames(P, x$nodes)
}

#' Pairwise layer similarity
#'
#' Cosine similarity between the (upper-triangular) edge-weight vectors of each
#' pair of layers: 1 = identical coupling structure, 0 = disjoint.
#'
#' @param x A `"multiplex"` object.
#' @return A `n_layers x n_layers` similarity matrix.
#' @export
layerSimilarity <- function(x) {
  stopifnot(inherits(x, "multiplex"))
  ut <- upper.tri(x$layers[[1]])
  vecs <- vapply(x$layers, function(m) m[ut], numeric(sum(ut)))
  S <- matrix(NA_real_, x$n_layers, x$n_layers,
              dimnames = list(x$layer_names, x$layer_names))
  for (a in seq_len(x$n_layers)) for (b in seq_len(x$n_layers)) {
    va <- vecs[, a]; vb <- vecs[, b]
    den <- sqrt(sum(va^2)) * sqrt(sum(vb^2))
    S[a, b] <- if (den > 0) sum(va * vb) / den else 0
  }
  S
}

#' Aggregate a multiplex into a single weighted network
#'
#' @param x A `"multiplex"` object.
#' @param method `"sum"` (default) or `"mean"` over layers.
#' @return A single weighted adjacency matrix.
#' @export
aggregateNetwork <- function(x, method = c("sum", "mean")) {
  stopifnot(inherits(x, "multiplex"))
  method <- match.arg(method)
  A <- Reduce(`+`, x$layers)
  if (method == "mean") A <- A / x$n_layers
  A
}

#' Supra-adjacency matrix of a multiplex
#'
#' Builds the `(N*L) x (N*L)` block matrix with the layer adjacencies on the
#' diagonal blocks and `interlayer` coupling on the identity off-diagonal blocks
#' (categorical coupling of a node to its own replicas across layers).
#'
#' @param x A `"multiplex"` object.
#' @param interlayer Interlayer coupling weight (default 1).
#' @return A `(N*L) x (N*L)` numeric matrix.
#' @export
supraAdjacency <- function(x, interlayer = 1) {
  stopifnot(inherits(x, "multiplex"))
  N <- length(x$nodes); L <- x$n_layers
  S <- matrix(0, N * L, N * L)
  for (a in seq_len(L)) {
    r <- ((a - 1) * N + 1):(a * N)
    S[r, r] <- x$layers[[a]]
  }
  if (interlayer != 0) for (a in seq_len(L)) for (b in seq_len(L)) if (a != b) {
    ra <- ((a - 1) * N + 1):(a * N)
    rb <- ((b - 1) * N + 1):(b * N)
    S[cbind(ra, rb)] <- interlayer
  }
  S
}

#' @export
print.multiplex <- function(x, ...) {
  cat("<Multiplex physiological network>\n")
  cat(sprintf("  nodes:  %d (%s)\n", length(x$nodes),
              paste(utils::head(x$nodes, 6), collapse = ", ")))
  cat(sprintf("  layers: %d (%s)\n", x$n_layers,
              paste(x$layer_names, collapse = ", ")))
  invisible(x)
}
