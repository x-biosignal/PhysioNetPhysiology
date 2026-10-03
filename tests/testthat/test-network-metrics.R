library(testthat)
library(PhysioNetPhysiology)

# a clear two-block network with a weak inter-block bridge
two_block <- function() {
  A <- matrix(0, 10, 10)
  A[1:5, 1:5] <- 0.9; A[6:10, 6:10] <- 0.9
  A[1, 6] <- A[6, 1] <- 0.15
  diag(A) <- 0
  rownames(A) <- colnames(A) <- paste0("n", 1:10)
  A
}

test_that("communityDetection recovers a two-block partition with high modularity", {
  skip_if_not_installed("igraph")
  cd <- communityDetection(two_block(), seed = 1)
  expect_s3_class(cd, "community_partition")
  expect_equal(cd$n_communities, 2L)
  expect_gt(cd$modularity, 0.3)
  # the two blocks are separated
  m <- cd$membership
  expect_equal(length(unique(m[1:5])), 1L)
  expect_equal(length(unique(m[6:10])), 1L)
  expect_false(m[["n1"]] == m[["n6"]])
})

test_that("communityDetection membership feeds the existing role metrics", {
  skip_if_not_installed("igraph")
  A <- two_block()
  cd <- communityDetection(A, seed = 1)
  roles <- moduleRoles(A, cd$membership)          # closes the loop role metrics needed
  expect_true(is.data.frame(roles) || is.list(roles))
})

test_that("networkMetrics returns sane global + node metrics (matrix and tds_network)", {
  skip_if_not_installed("igraph")
  A <- two_block()
  nm <- networkMetrics(A)
  expect_s3_class(nm, "network_metrics")
  g <- nm$global
  expect_true(g$density >= 0 && g$density <= 1)
  expect_true(g$transitivity >= 0 && g$transitivity <= 1)
  expect_true(g$global_efficiency >= 0 && g$global_efficiency <= 1)
  expect_equal(g$n_communities, 2L)
  expect_setequal(names(nm$node),
    c("node", "degree", "strength", "clustering", "betweenness",
      "closeness", "eigenvector", "local_efficiency"))
  expect_true(all(nm$node$clustering >= 0 & nm$node$clustering <= 1))
  expect_true(all(nm$node$local_efficiency >= 0 & nm$node$local_efficiency <= 1))

  # accepts a tds_network object (uses its $adjacency)
  net <- structure(list(adjacency = A), class = "tds_network")
  nm2 <- networkMetrics(net)
  expect_equal(nm2$global$n_nodes, 10L)
})

test_that("smallWorldness flags a Watts-Strogatz graph and not a random one", {
  skip_if_not_installed("igraph")
  gws <- igraph::sample_smallworld(1, 60, 4, 0.08)
  grn <- igraph::sample_gnm(60, igraph::gsize(gws))
  sw <- smallWorldness(as.matrix(igraph::as_adjacency_matrix(gws)), n_rand = 25, seed = 2)
  sr <- smallWorldness(as.matrix(igraph::as_adjacency_matrix(grn)), n_rand = 25, seed = 2)
  expect_s3_class(sw, "small_world")
  expect_gt(sw$sigma, 1.5)                         # small-world
  expect_gt(sw$sigma, sr$sigma)                    # more small-world than random
})

test_that("network functions validate input and require igraph", {
  expect_error(PhysioNetPhysiology:::.pn_adjacency(matrix(0, 2, 3)), "square")
  if (requireNamespace("igraph", quietly = TRUE)) {
    expect_error(networkMetrics("not a network"), "adjacency matrix or a tds_network")
  }
})
