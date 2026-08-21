# Graph-topology metrics and community detection for physiological networks.
#
# Classical network-science measures (clustering, path length, efficiency,
# centralities, assortativity, modularity, community structure, small-worldness)
# on the weighted networks produced by tdsNetwork()/tdsNetworkByState(). These
# delegate to the validated `igraph` library (a package Suggests) -- the graph
# algorithms (Dijkstra, betweenness, Louvain, ...) are its remit, not something
# to re-derive. communityDetection() also supplies the partition that the
# existing role metrics (participationCoefficient / moduleRoles) require.

.need_igraph <- function(fn) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("The 'igraph' package is required for ", fn,
         "(). Install it with install.packages('igraph').", call. = FALSE)
}

# adjacency matrix (zero diagonal) from a tds_network or a plain matrix
.pn_adjacency <- function(x) {
  adj <- if (inherits(x, "tds_network")) x$adjacency
         else if (is.matrix(x) || is.data.frame(x)) as.matrix(x)
         else stop("x must be a numeric adjacency matrix or a tds_network.", call. = FALSE)
  if (nrow(adj) != ncol(adj)) stop("adjacency must be square.", call. = FALSE)
  storage.mode(adj) <- "double"
  diag(adj) <- 0
  adj
}

.pn_node_names <- function(adj) {
  rn <- rownames(adj)
  if (is.null(rn)) paste0("V", seq_len(nrow(adj))) else rn
}

.pn_igraph <- function(adj, directed = FALSE) {
  igraph::graph_from_adjacency_matrix(
    adj, mode = if (directed) "directed" else "undirected",
    weighted = TRUE, diag = FALSE)
}

#' Detect communities in a physiological network
#'
#' Partitions a network into communities (modules) and reports the modularity Q.
#' The membership can be passed to [moduleRoles()] / [participationCoefficient()],
#' closing the loop those role metrics need. Delegates to `igraph`.
#'
#' @param x A `tds_network` (from [tdsNetwork()]) or a weighted adjacency matrix.
#' @param method Community-detection algorithm: `"louvain"` (default),
#'   `"fast_greedy"`, `"walktrap"`, `"leading_eigen"`, or `"infomap"`.
#' @param seed Optional integer seed (some algorithms are stochastic).
#' @return A `community_partition` object: `membership` (named integer vector),
#'   `modularity` (Q), `n_communities`, `sizes`, and `method`.
#' @seealso [networkMetrics()], [moduleRoles()]
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' A <- matrix(0, 8, 8)
#' A[1:4, 1:4] <- 0.9; A[5:8, 5:8] <- 0.9; A[1, 5] <- A[5, 1] <- 0.2
#' diag(A) <- 0
#' communityDetection(A)
#' }
communityDetection <- function(x, method = c("louvain", "fast_greedy",
                                             "walktrap", "leading_eigen", "infomap"),
                               seed = NULL) {
  .need_igraph("communityDetection")
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  adj <- .pn_adjacency(x)
  nodes <- .pn_node_names(adj)
  g <- .pn_igraph(adj, directed = FALSE)          # modularity communities are undirected
  w <- igraph::E(g)$weight
  comm <- switch(method,
    louvain       = igraph::cluster_louvain(g, weights = w),
    fast_greedy   = igraph::cluster_fast_greedy(g, weights = w),
    walktrap      = igraph::cluster_walktrap(g, weights = w),
    leading_eigen = igraph::cluster_leading_eigen(g, weights = w),
    infomap       = igraph::cluster_infomap(g, e.weights = w))
  mem <- stats::setNames(as.integer(igraph::membership(comm)), nodes)
  structure(list(membership = mem, modularity = igraph::modularity(comm),
                 n_communities = length(unique(mem)),
                 sizes = as.integer(igraph::sizes(comm)), method = method),
            class = "community_partition")
}

#' @export
print.community_partition <- function(x, ...) {
  cat(sprintf("<community_partition> method = %s\n", x$method))
  cat(sprintf("  %d communities (sizes %s); modularity Q = %.3f\n",
              x$n_communities, paste(x$sizes, collapse = ", "), x$modularity))
  invisible(x)
}

#' Graph-topology metrics of a physiological network
#'
#' Computes the standard network-science metrics on a weighted physiological
#' network: node-level (degree, strength, weighted clustering, betweenness,
#' closeness, eigenvector centrality, local efficiency) and global (density,
#' transitivity, characteristic path length, global/local efficiency, degree
#' assortativity, and modularity from a community partition). Edge weights are
#' treated as coupling strengths; path-based metrics use distances `1/weight`
#' (stronger coupling = shorter path). Delegates to `igraph`.
#'
#' @param x A `tds_network` or a weighted adjacency matrix.
#' @param directed Treat the network as directed (default `FALSE`).
#' @param community_method Community algorithm for the modularity summary
#'   (see [communityDetection()]).
#' @return A `network_metrics` object: `global` (named list), `node` (data frame
#'   of per-node metrics), `membership`, and `directed`.
#' @seealso [communityDetection()], [smallWorldness()], [moduleRoles()]
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' A <- matrix(runif(64, 0, 0.3), 8, 8); A <- (A + t(A)) / 2; diag(A) <- 0
#' nm <- networkMetrics(A)
#' nm$global$global_efficiency
#' }
networkMetrics <- function(x, directed = FALSE, community_method = "louvain") {
  .need_igraph("networkMetrics")
  adj <- .pn_adjacency(x)
  n <- nrow(adj)
  nodes <- .pn_node_names(adj)
  g <- .pn_igraph(adj, directed)
  w <- igraph::E(g)$weight
  dist <- 1 / w                                   # coupling strength -> distance

  deg  <- igraph::degree(g)
  strg <- igraph::strength(g, weights = w)
  clus <- igraph::transitivity(g, type = "weighted", weights = w, isolates = "zero")
  clus <- pmin(pmax(clus, 0), 1)                  # Barrat clustering is bounded [0,1]
  betw <- igraph::betweenness(g, weights = dist, directed = directed)
  clos <- suppressWarnings(igraph::closeness(g, weights = dist))
  eig  <- igraph::eigen_centrality(g, weights = w, directed = directed)$vector
  leff <- igraph::local_efficiency(g, weights = dist, directed = directed)

  node_df <- data.frame(
    node = nodes, degree = deg, strength = strg, clustering = clus,
    betweenness = betw, closeness = clos, eigenvector = eig,
    local_efficiency = leff, row.names = NULL, stringsAsFactors = FALSE)

  comm <- communityDetection(adj, method = community_method)
  global <- list(
    n_nodes = n, n_edges = igraph::gsize(g),
    density = igraph::edge_density(g),
    mean_degree = mean(deg), mean_strength = mean(strg),
    transitivity = igraph::transitivity(g, type = "global"),
    char_path_length = igraph::mean_distance(g, weights = dist, directed = directed),
    global_efficiency = igraph::global_efficiency(g, weights = dist, directed = directed),
    mean_local_efficiency = mean(leff, na.rm = TRUE),
    assortativity = tryCatch(igraph::assortativity_degree(g, directed = directed),
                             error = function(e) NA_real_),
    modularity = comm$modularity, n_communities = comm$n_communities)

  structure(list(global = global, node = node_df,
                 membership = comm$membership, directed = directed),
            class = "network_metrics")
}

#' @export
print.network_metrics <- function(x, ...) {
  gl <- x$global
  cat(sprintf("<network_metrics> %d nodes, %d edges%s\n", gl$n_nodes, gl$n_edges,
              if (x$directed) " (directed)" else ""))
  cat(sprintf("  density %.3f | transitivity %.3f | char path %.3f | global eff %.3f\n",
              gl$density, gl$transitivity, gl$char_path_length, gl$global_efficiency))
  cat(sprintf("  assortativity %.3f | modularity %.3f (%d communities)\n",
              gl$assortativity, gl$modularity, gl$n_communities))
  invisible(x)
}

#' Small-worldness of a physiological network
#'
#' Quantifies small-world organization on the binarized network: the
#' Humphries-Gurney sigma = (C/C_rand)/(L/L_rand) (small-world when > 1, comparing
#' clustering and path length to degree-preserving random graphs) and the Telford
#' omega = L_rand/L - C/C_latt (near 0 = small-world; -> 1 random; -> -1 lattice).
#' Delegates the randomization and lattice references to `igraph`.
#'
#' @param x A `tds_network` or a weighted adjacency matrix (binarized by
#'   presence of an edge).
#' @param n_rand Number of degree-preserving random graphs (default 50).
#' @param seed Optional integer seed for the randomization.
#' @return A `small_world` object: `sigma`, `omega`, and the components `C`, `L`,
#'   `C_rand`, `L_rand`, `C_latt`, `n_rand`.
#' @references Humphries & Gurney 2008; Telford et al. 2011.
#' @seealso [networkMetrics()]
#' @export
#' @examples
#' \donttest{
#' g <- igraph::sample_smallworld(1, 40, 4, 0.05)     # a small-world graph
#' sw <- smallWorldness(as.matrix(igraph::as_adjacency_matrix(g)), n_rand = 20, seed = 1)
#' sw$sigma
#' }
smallWorldness <- function(x, n_rand = 50L, seed = NULL) {
  .need_igraph("smallWorldness")
  if (!is.null(seed)) set.seed(seed)
  adj <- .pn_adjacency(x)
  gb <- igraph::graph_from_adjacency_matrix(adj > 0, mode = "undirected", diag = FALSE)
  C <- igraph::transitivity(gb, type = "global")
  L <- igraph::mean_distance(gb)
  deg <- igraph::degree(gb)

  Cr <- numeric(n_rand); Lr <- numeric(n_rand)
  for (i in seq_len(n_rand)) {
    gr <- tryCatch(igraph::sample_degseq(deg, method = "vl"),
                   error = function(e)
                     igraph::rewire(gb, igraph::keeping_degseq(niter = 10 * igraph::gsize(gb))))
    Cr[i] <- igraph::transitivity(gr, type = "global")
    Lr[i] <- igraph::mean_distance(gr)
  }
  Crand <- mean(Cr, na.rm = TRUE); Lrand <- mean(Lr, na.rm = TRUE)
  sigma <- (C / Crand) / (L / Lrand)

  k <- max(1L, round(mean(deg)) %/% 2L)
  Clatt <- tryCatch(igraph::transitivity(
    igraph::sample_smallworld(1, nrow(adj), k, 0), type = "global"),
    error = function(e) NA_real_)
  omega <- if (is.finite(Clatt) && Clatt > 0) Lrand / L - C / Clatt else NA_real_

  structure(list(sigma = sigma, omega = omega, C = C, L = L,
                 C_rand = Crand, L_rand = Lrand, C_latt = Clatt, n_rand = n_rand),
            class = "small_world")
}

#' @export
print.small_world <- function(x, ...) {
  cat(sprintf("<small_world> sigma = %.3f  omega = %.3f\n", x$sigma, x$omega))
  cat(sprintf("  C = %.3f (rand %.3f, latt %.3f) | L = %.3f (rand %.3f)\n",
              x$C, x$C_rand, x$C_latt, x$L, x$L_rand))
  cat(sprintf("  %s\n", if (isTRUE(x$sigma > 1)) "small-world (sigma > 1)" else
    "not small-world (sigma <= 1)"))
  invisible(x)
}
