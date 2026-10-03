test_that("TDS detects a fixed-delay coupling and recovers its lag", {
  X <- make_lagged(delay = 5)
  res <- timeDelayStability(X, window = 100, step = 50, max_lag = 20,
                            min_stable = 3)

  expect_s3_class(res, "tds")
  expect_equal(dim(res$tds), c(3, 3))
  expect_equal(diag(res$tds), c(x = 0, y = 0, z = 0))

  # coupled pair is strongly stable; unrelated pairs are weak
  expect_gt(res$tds["x", "y"], 50)
  expect_lt(res$tds["x", "z"], 25)
  expect_lt(res$tds["y", "z"], 25)
  expect_gt(res$tds["x", "y"], res$tds["x", "z"] + 20)

  # recovered stable lag matches the planted delay (x leads y by 5)
  expect_equal(abs(res$stable_lag["x", "y"]), 5, tolerance = 2)
})

test_that("larger min_stable never increases TDS", {
  X <- make_lagged()
  t3 <- timeDelayStability(X, window = 100, step = 50, max_lag = 20,
                           min_stable = 3)$tds["x", "y"]
  t8 <- timeDelayStability(X, window = 100, step = 50, max_lag = 20,
                           min_stable = 8)$tds["x", "y"]
  expect_lte(t8, t3)
})

test_that("stronger noise weakens a coupling", {
  clean <- timeDelayStability(make_lagged(sd_noise = 0.2), window = 100,
                              step = 50, max_lag = 20, min_stable = 3)$tds["x", "y"]
  noisy <- timeDelayStability(make_lagged(sd_noise = 2.0), window = 100,
                              step = 50, max_lag = 20, min_stable = 3)$tds["x", "y"]
  expect_gt(clean, noisy)
})

test_that("summary() ranks pairs by TDS", {
  res <- timeDelayStability(make_lagged(), window = 100, step = 50,
                            max_lag = 20, min_stable = 3)
  s <- summary(res)
  expect_s3_class(s, "data.frame")
  expect_equal(s$from[1], "x")
  expect_equal(s$to[1], "y")             # strongest pair first
  expect_true(all(diff(s$tds) <= 0))     # descending
})

test_that("summary() keeps pair labels aligned with values at >=4 nodes", {
  # regression: upper.tri (column-major) vs combn ordering diverge for n >= 4.
  # Build 4 nodes where only node 2 (b) and node 3 (c) are coupled.
  set.seed(11)
  n <- 3000
  a <- rnorm(n)
  b <- as.numeric(stats::filter(rnorm(n), rep(1 / 5, 5), sides = 2)); b[is.na(b)] <- 0
  cc <- c(rep(0, 5), b[seq_len(n - 5)]) + rnorm(n, sd = 0.3)   # c lags b by 5
  d <- rnorm(n)
  X <- cbind(a = a, b = b, c = cc, d = d)
  res <- timeDelayStability(X, window = 100, step = 50, max_lag = 20,
                            min_stable = 3)
  s <- summary(res)
  # every summary row must equal the matrix entry it names
  for (i in seq_len(nrow(s)))
    expect_equal(s$tds[i], res$tds[s$from[i], s$to[i]])
  # the coupled pair is b-c, and it is the strongest
  expect_setequal(c(s$from[1], s$to[1]), c("b", "c"))
})

test_that("input validation", {
  X <- make_lagged(n = 200)
  expect_error(timeDelayStability(X, window = 500), "longer than")
  expect_error(timeDelayStability(X[, 1, drop = FALSE], window = 50),
               "at least two")
  Xna <- X; Xna[1, 1] <- NA
  expect_error(timeDelayStability(Xna, window = 50), "NA")
  expect_error(timeDelayStability(X, window = 50, max_lag = 60), "max_lag")
})
