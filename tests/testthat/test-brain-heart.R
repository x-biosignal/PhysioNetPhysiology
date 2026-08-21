make_hep_data <- function(n_sec = 60, fs = 250, amp = 3, seed = 1) {
  set.seed(seed)
  n <- n_sec * fs
  rpeaks <- seq(fs, n - fs, by = round(0.9 * fs))     # ~0.9 s beats
  eeg <- stats::rnorm(n)
  # a Gaussian deflection 0.3 s after each R-peak
  bt <- round(0.3 * fs); w <- round(0.04 * fs)
  for (r in rpeaks) {
    idx <- (r + bt - 3 * w):(r + bt + 3 * w)
    idx <- idx[idx >= 1 & idx <= n]
    eeg[idx] <- eeg[idx] + amp * exp(-((idx - (r + bt))^2) / (2 * w^2))
  }
  list(eeg = eeg, rpeaks = rpeaks, fs = fs)
}

test_that("HEP recovers a heartbeat-locked deflection", {
  d <- make_hep_data()
  hep <- heartbeatEvokedPotential(d$eeg, d$rpeaks, d$fs, window = c(-0.1, 0.6))

  expect_s3_class(hep, "hep")
  expect_equal(hep$n, length(d$rpeaks) - 1, tolerance = 1)  # edge epochs dropped
  pk <- hep$time[which.max(hep$hep)]
  expect_equal(pk, 0.3, tolerance = 0.05)                   # peak near 300 ms
  expect_gt(max(hep$hep), 1)                                # visible after averaging
})

test_that("HEP amplitude reader targets the component window", {
  d <- make_hep_data()
  hep <- heartbeatEvokedPotential(d$eeg, d$rpeaks, d$fs)
  a <- hepAmplitude(hep, component = c(0.25, 0.35), fun = "peak")
  expect_gt(a, 1)
  # a window with no locked activity is near zero
  b <- hepAmplitude(hep, component = c(-0.1, -0.05), fun = "mean")
  expect_lt(abs(b), 0.5)
})

test_that("permutation test flags the real HEP and not shuffled peaks", {
  d <- make_hep_data(amp = 3)
  sig <- hepSignificance(d$eeg, d$rpeaks, d$fs, component = c(0.25, 0.35),
                         n_perm = 200)
  expect_lt(sig$p_value, 0.05)
  expect_gt(sig$observed, stats::quantile(sig$null, 0.95))

  # a brain signal with no cardiac locking should not be significant
  set.seed(9)
  eeg0 <- stats::rnorm(length(d$eeg))
  sig0 <- hepSignificance(eeg0, d$rpeaks, d$fs, component = c(0.25, 0.35),
                          n_perm = 200)
  expect_gt(sig0$p_value, 0.05)
})

test_that("HEP errors when the window has too few usable beats", {
  expect_error(
    heartbeatEvokedPotential(rnorm(100), c(5), 250, window = c(-0.1, 0.6)),
    "too few")
})
