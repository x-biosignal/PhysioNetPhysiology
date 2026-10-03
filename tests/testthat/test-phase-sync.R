test_that("phase sync index: identical = 1, constant-shift ~ 1, independent low", {
  fs <- 100; t <- seq(0, 60, by = 1 / fs)
  pa  <- instantaneousPhase(sin(2 * pi * 0.5 * t), fs)
  pa2 <- instantaneousPhase(sin(2 * pi * 0.5 * t + 0.7), fs)   # constant offset
  pb  <- instantaneousPhase(sin(2 * pi * 0.37 * t), fs)        # different freq

  expect_equal(phaseSynchronizationIndex(pa, pa), 1, tolerance = 1e-6)
  expect_gt(phaseSynchronizationIndex(pa, pa2), 0.95)
  expect_lt(phaseSynchronizationIndex(pa, pb), 0.6)
})

test_that("n:m scan recovers a 4:1 frequency ratio", {
  fs <- 100; t <- seq(0, 120, by = 1 / fs)
  pr <- instantaneousPhase(sin(2 * pi * 0.25 * t), fs)   # resp 0.25 Hz
  pc <- instantaneousPhase(sin(2 * pi * 1.00 * t), fs)   # cardiac 1 Hz (4x)
  nm <- nmPhaseSynchronization(pc, pr, orders_n = 1:2, orders_m = 1:8)

  expect_equal(nm$ratio, "4:1")          # f_cardiac : f_resp
  expect_gt(nm$lambda, 0.9)
})

test_that("cardiacPhaseFromPeaks resets each beat", {
  rp <- seq(10, 500, by = 25)
  ph <- cardiacPhaseFromPeaks(rp, 600)
  vals <- ph[!is.na(ph)]
  expect_true(all(vals >= 0 & vals < 2 * pi))
  expect_true(is.na(ph[1]))              # before the first peak
  expect_lt(ph[rp[2]], 0.3)              # ~0 at a beat mark
})

test_that("cardiorespiratoryCoupling finds 4 heartbeats per breath", {
  fs <- 100; t <- seq(0, 120, by = 1 / fs)
  resp <- sin(2 * pi * 0.25 * t)
  card <- sin(2 * pi * 1.00 * t)
  rpeaks <- which(diff(sign(card)) > 0)  # upward zero crossings ~ 1 Hz

  cr <- cardiorespiratoryCoupling(resp, rpeaks, fs)
  expect_s3_class(cr, "cardioresp")
  expect_equal(cr$ratio, "4:1")
  expect_gt(cr$lambda, 0.8)
  expect_s3_class(cr$synchrogram, "synchrogram")
  expect_equal(nrow(cr$synchrogram), length(rpeaks))
  expect_true(all(cr$synchrogram$psi >= 0 & cr$synchrogram$psi < 1))
})

test_that("bandpass isolates a component before phase extraction", {
  fs <- 100; t <- seq(0, 30, by = 1 / fs)
  x <- sin(2 * pi * 0.3 * t) + 2 * sin(2 * pi * 10 * t)   # slow + fast
  # unbandpassed phase is dominated by the 10 Hz term; banded recovers 0.3 Hz
  ph <- instantaneousPhase(x, fs, band = c(0.1, 0.6))
  ref <- instantaneousPhase(sin(2 * pi * 0.3 * t), fs)
  # compare on the interior (avoid FFT edge effects)
  keep <- seq(300, length(t) - 300)
  expect_gt(phaseSynchronizationIndex(ph[keep], ref[keep]), 0.9)
})
