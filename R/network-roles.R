# Module-role metrics (Guimera & Amaral 2005) -----------------------------
# Given a network and a community partition, characterise each node's role by
# how its links distribute within vs across modules: the participation
# coefficient (across-module spread) and the within-module degree z-score.

#' @keywords internal
#' @noRd
.check_roles_input <- function(adj, communities) {
  adj <- as.matrix(adj)
  if (nrow(adj) != ncol(adj))
    stop("`adj` must be square.", call. = FALSE)
  if (length(communities) != nrow(adj))
    stop("`communities` must have one entry per node.", call. = FALSE)
  invisible(adj)
}

#' Participation coefficient
#'
#' \eqn{P_i = 1 - \sum_s (\kappa_{is}/k_i)^2}, where \eqn{\kappa_{is}} is node
#' `i`'s summed connection weight to module `s` and \eqn{k_i} its total strength.
#' 0 = all links inside one module; approaches 1 as links spread evenly over
#' many modules.
#'
#' @param adj A square (weighted or binary) adjacency matrix.
#' @param communities A vector of module labels, one per node.
#' @return A named numeric vector of participation coefficients in `[0, 1)`.
#' @examples
#' adj <- matrix(0, 4, 4); adj[1,2] <- adj[2,1] <- 1; adj[3,4] <- adj[4,3] <- 1
#' participationCoefficient(adj, c(1, 1, 2, 2))
#' @export
participationCoefficient <- function(adj, communities) {
  adj <- .check_roles_input(adj, communities)
  ki <- rowSums(adj)
  mods <- unique(communities)
  ksum <- vapply(mods, function(s)
    rowSums(adj[, communities == s, drop = FALSE]), numeric(nrow(adj)))
  P <- 1 - rowSums((ksum / ki)^2)
  P[ki == 0] <- 0
  stats::setNames(P, rownames(adj) %||% seq_len(nrow(adj)))
}

#' Within-module degree z-score
#'
#' Standardises each node's within-module strength relative to the other nodes
#' in its module. High values mark within-module hubs.
#'
#' @param adj A square adjacency matrix.
#' @param communities A vector of module labels, one per node.
#' @return A named numeric vector of z-scores.
#' @export
withinModuleDegreeZ <- function(adj, communities) {
  adj <- .check_roles_input(adj, communities)
  z <- numeric(nrow(adj))
  for (s in unique(communities)) {
    idx <- which(communities == s)
    kin <- rowSums(adj[idx, idx, drop = FALSE])
    mu <- mean(kin); sdv <- stats::sd(kin)
    z[idx] <- if (is.finite(sdv) && sdv > 0) (kin - mu) / sdv else 0
  }
  stats::setNames(z, rownames(adj) %||% seq_len(nrow(adj)))
}

#' Node module roles (Guimera-Amaral)
#'
#' Combines the within-module degree z-score and participation coefficient into
#' the seven-role classification of Guimera & Amaral (2005).
#'
#' @param adj A square adjacency matrix.
#' @param communities A vector of module labels, one per node.
#' @return A data frame with `node`, `z` (within-module degree z-score), `P`
#'   (participation coefficient), and `role`.
#' @export
moduleRoles <- function(adj, communities) {
  z <- withinModuleDegreeZ(adj, communities)
  P <- participationCoefficient(adj, communities)
  role <- mapply(function(zi, pi) {
    if (zi >= 2.5) {                              # hubs
      if (pi < 0.30) "provincial hub"
      else if (pi < 0.75) "connector hub"
      else "kinless hub"
    } else {                                      # non-hubs
      if (pi < 0.05) "ultra-peripheral"
      else if (pi < 0.62) "peripheral"
      else if (pi < 0.80) "connector"
      else "kinless"
    }
  }, z, P)
  data.frame(node = names(z), z = z, P = P, role = role, row.names = NULL)
}
