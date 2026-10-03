# Shared synthetic-data generators for the TDS test suite.

# A smoothed (autocorrelated) driver x, a copy y that lags x by `delay`
# samples plus noise, and an unrelated white-noise node z.
make_lagged <- function(n = 3000, delay = 5, sd_noise = 0.3, seed = 1) {
  set.seed(seed)
  x <- as.numeric(stats::filter(stats::rnorm(n), rep(1 / 5, 5), sides = 2))
  x[is.na(x)] <- 0
  y <- c(rep(0, delay), x[seq_len(n - delay)]) + stats::rnorm(n, sd = sd_noise)
  z <- stats::rnorm(n)
  cbind(x = x, y = y, z = z)
}

# Two-state series: node `heart` couples to `brain` (delay 5) only in state B.
make_two_state <- function(n = 4000, delay = 5, sd_noise = 0.3, seed = 5) {
  set.seed(seed)
  half <- n %/% 2
  x <- as.numeric(stats::filter(stats::rnorm(n), rep(1 / 5, 5), sides = 2))
  x[is.na(x)] <- 0
  y <- stats::rnorm(n, sd = sd_noise)
  # y(t) = x(t - delay) for t in the second half (state B)
  y[(half + delay + 1):n] <- x[(half + 1):(n - delay)] +
    stats::rnorm(n - half - delay, sd = sd_noise)
  list(X = cbind(brain = x, heart = y),
       states = rep(c("A", "B"), each = half))
}
