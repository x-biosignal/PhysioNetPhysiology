test_that("C++ peak-lag engine matches the pure-R reference", {
  set.seed(1)
  n <- 3000
  x <- rnorm(n)
  y <- c(rep(0, 7), x[seq_len(n - 7)]) + rnorm(n, sd = 0.4)
  win <- 200L; step <- 50L; max_lag <- 30L
  starts <- seq.int(1L, n - win + 1L, by = step)

  cpp <- tds_peaklag_cpp(as.double(x), as.double(y),
                         as.integer(starts), win, max_lag)
  ref <- .pair_peaklag_r(x, y, starts, win, max_lag)

  expect_identical(cpp$tau0, as.integer(ref$tau0))     # identical optimal lags
  expect_equal(cpp$peak, ref$peak, tolerance = 1e-9)   # identical peak cc
})

test_that("timeDelayStability still detects a known coupling (C++ path)", {
  X <- make_lagged(delay = 5)
  res <- timeDelayStability(X, window = 100, step = 50, max_lag = 20,
                            min_stable = 3)
  expect_gt(res$tds["x", "y"], 50)
  expect_lt(res$tds["x", "z"], 25)
  expect_equal(abs(res$stable_lag["x", "y"]), 5, tolerance = 2)
})
