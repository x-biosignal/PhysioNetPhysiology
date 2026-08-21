test_that("network thresholding keeps strong links and drops weak ones", {
  res <- timeDelayStability(make_lagged(), window = 100, step = 50,
                            max_lag = 20, min_stable = 3)
  net <- tdsNetwork(res, threshold = 40)

  expect_s3_class(net, "tds_network")
  expect_true(any(net$edges$from == "x" & net$edges$to == "y"))
  expect_equal(unname(net$adjacency["x", "z"]), 0)
  expect_equal(unname(net$degree[["x"]]), 1)
  expect_equal(net$n_links, 1L)
})

test_that("edges carry the stable lag", {
  res <- timeDelayStability(make_lagged(), window = 100, step = 50,
                            max_lag = 20, min_stable = 3)
  net <- tdsNetwork(res, threshold = 40)
  expect_equal(abs(net$edges$lag[1]), 5, tolerance = 2)
})

test_that("surrogate threshold sits below the real coupling", {
  set.seed(3)
  X <- make_lagged()
  thr <- tdsSurrogateThreshold(X, n_surrogates = 20, window = 100, step = 50,
                               max_lag = 20, min_stable = 3)
  expect_true(thr > 0 && thr < 100)

  res <- timeDelayStability(X, window = 100, step = 50, max_lag = 20,
                            min_stable = 3)
  expect_gt(res$tds["x", "y"], thr)          # real coupling clears the floor
  expect_lt(thr, res$tds["x", "y"])
})

test_that("tdsNetwork can derive its own surrogate threshold", {
  set.seed(4)
  X <- make_lagged()
  res <- timeDelayStability(X, window = 100, step = 50, max_lag = 20,
                            min_stable = 3)
  net <- tdsNetwork(res, surrogate = TRUE, X = X, n_surrogates = 15)
  expect_s3_class(net, "tds_network")
  expect_gt(net$adjacency["x", "y"], 0)
})

test_that("tdsNetwork errors without a threshold or surrogate", {
  res <- timeDelayStability(make_lagged(), window = 100, step = 50,
                            max_lag = 20, min_stable = 3)
  expect_error(tdsNetwork(res), "threshold")
})
