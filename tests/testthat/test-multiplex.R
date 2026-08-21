make_layers <- function() {
  nodes <- c("brain", "heart", "resp", "muscle")
  L1 <- matrix(0, 4, 4, dimnames = list(nodes, nodes))
  L1["brain", "heart"] <- L1["heart", "brain"] <- 1
  L1["heart", "resp"]  <- L1["resp", "heart"]  <- 1
  L2 <- matrix(0, 4, 4, dimnames = list(nodes, nodes))
  L2["brain", "heart"] <- L2["heart", "brain"] <- 1     # shared with L1
  L2["brain", "muscle"] <- L2["muscle", "brain"] <- 1   # unique to L2
  list(alpha = L1, beta = L2)
}

test_that("multiplexNetwork validates and stores layers", {
  mp <- multiplexNetwork(make_layers())
  expect_s3_class(mp, "multiplex")
  expect_equal(mp$n_layers, 2)
  expect_equal(mp$nodes, c("brain", "heart", "resp", "muscle"))
  expect_error(multiplexNetwork(list(matrix(0, 2, 2))), "at least two")
  expect_error(
    multiplexNetwork(list(a = matrix(0, 3, 3), b = matrix(0, 2, 2))),
    "same dimensions")
})

test_that("overlapping degree sums per-layer degree", {
  d <- multiplexDegree(multiplexNetwork(make_layers()))
  # brain: L1 -> heart (deg1); L2 -> heart, muscle (deg2) => overlapping 3
  expect_equal(unname(d$overlapping["brain"]), 3)
  expect_equal(unname(d$per_layer["brain", "beta"]), 2)
})

test_that("participation is 0 for a single-layer node, higher when spread", {
  mp <- multiplexNetwork(make_layers())
  P <- multiplexParticipation(mp)
  # resp only appears in layer alpha -> participation 0
  expect_equal(unname(P["resp"]), 0)
  # heart spreads across both layers -> participation > 0
  expect_gt(P["heart"], 0)
  expect_true(all(P >= 0 & P <= 1))
})

test_that("layer similarity is 1 on the diagonal and <1 for distinct layers", {
  S <- layerSimilarity(multiplexNetwork(make_layers()))
  expect_equal(diag(S), c(alpha = 1, beta = 1), tolerance = 1e-8)
  expect_lt(S["alpha", "beta"], 1)
  expect_gt(S["alpha", "beta"], 0)          # they share the brain-heart edge
})

test_that("aggregate and supra-adjacency have the right shape", {
  mp <- multiplexNetwork(make_layers())
  A <- aggregateNetwork(mp, "sum")
  expect_equal(unname(A["brain", "heart"]), 2)   # edge present in both layers
  S <- supraAdjacency(mp, interlayer = 1)
  expect_equal(dim(S), c(8, 8))                  # 4 nodes x 2 layers
  expect_equal(S[1, 5], 1)                       # brain(layer1) - brain(layer2)
})

test_that("module roles: participation and within-module z", {
  # two clean modules (1,2)-(3,4); node 2 also links to module 2 (a connector)
  adj <- matrix(0, 4, 4)
  adj[1, 2] <- adj[2, 1] <- 1
  adj[3, 4] <- adj[4, 3] <- 1
  adj[2, 3] <- adj[3, 2] <- 1                     # bridge between modules
  comm <- c(1, 1, 2, 2)

  P <- participationCoefficient(adj, comm)
  expect_equal(unname(P[1]), 0)                   # node 1 only links within
  expect_gt(P[2], 0)                              # node 2 bridges modules
  roles <- moduleRoles(adj, comm)
  expect_named(roles, c("node", "z", "P", "role"))
  expect_equal(nrow(roles), 4)
})
