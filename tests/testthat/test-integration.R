make_pe <- function(nt = 1000, fs = 100, chans = c("brain", "heart", "resp")) {
  set.seed(1)
  M <- vapply(seq_along(chans), function(k) {
    v <- as.numeric(stats::filter(stats::rnorm(nt), rep(1 / 5, 5), sides = 2))
    v[is.na(v)] <- 0; v
  }, numeric(nt))
  pe <- PhysioCore::PhysioExperiment(
    assays = list(raw = M), samplingRate = fs,
    colData = S4Vectors::DataFrame(name = chans))
  PhysioCore::channelNames(pe) <- chans
  pe
}

test_that("physioNodeMatrix extracts channels from a PhysioExperiment", {
  skip_if_not_installed("PhysioCore")
  pe <- make_pe()
  M <- physioNodeMatrix(pe)
  expect_true(is.matrix(M))
  expect_equal(ncol(M), 3)
  expect_equal(colnames(M), c("brain", "heart", "resp"))
  expect_equal(nrow(M), 1000)
})

test_that("envelope and down-sampling options work", {
  skip_if_not_installed("PhysioCore")
  pe <- make_pe()
  Me <- physioNodeMatrix(pe, feature = "envelope")
  expect_true(all(Me >= 0))
  expect_equal(colnames(Me), c("brain", "heart", "resp"))

  Md <- physioNodeMatrix(pe, target_rate = 10)     # 100 Hz -> 10 Hz
  expect_equal(nrow(Md), 100)
  expect_equal(ncol(Md), 3)
})

test_that("a PhysioExperiment flows straight into timeDelayStability", {
  skip_if_not_installed("PhysioCore")
  pe <- make_pe()
  res <- timeDelayStability(pe, window = 100, step = 50, max_lag = 20,
                            min_stable = 3)
  expect_s3_class(res, "tds")
  expect_equal(res$nodes, c("brain", "heart", "resp"))
})

test_that("physioNodeMatrix assembles nodes from a MultiPhysioExperiment", {
  skip_if_not_installed("PhysioCore")
  skip_if_not_installed("PhysioCrossModal")
  pe1 <- make_pe(nt = 1000, fs = 100, chans = c("eegband"))
  pe2 <- make_pe(nt = 500,  fs = 50,  chans = c("hr"))
  mpe <- PhysioCrossModal::MultiPhysioExperiment(
    experiments = list(eeg = pe1, ecg = pe2))
  M <- physioNodeMatrix(mpe, target_rate = 50)
  expect_true(is.matrix(M))
  expect_true(all(c("eeg.eegband", "ecg.hr") %in% colnames(M)))
  # differing rates without a target_rate is an error
  expect_error(physioNodeMatrix(mpe), "different sampling rates")
})
