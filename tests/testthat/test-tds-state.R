test_that("the organ-interaction network reconfigures across states", {
  d <- make_two_state()
  res <- tdsNetworkByState(d$X, d$states, window = 100, step = 50,
                           max_lag = 20, min_stable = 3, threshold = 30)

  expect_s3_class(res, "tds_state")
  expect_setequal(res$states, c("A", "B"))

  # coupling exists in B, not in A
  expect_gt(res$tds[["B"]]["brain", "heart"], res$tds[["A"]]["brain", "heart"])
  expect_equal(res$networks[["A"]]$n_links, 0L)
  expect_equal(res$networks[["B"]]$n_links, 1L)

  # reconfiguration summary reflects that
  recon <- res$reconfiguration
  expect_equal(recon$n_links[recon$state == "A"], 0L)
  expect_equal(recon$n_links[recon$state == "B"], 1L)
})

test_that("tdsReconfiguration reports a positive delta for the emerging link", {
  d <- make_two_state()
  res <- tdsNetworkByState(d$X, d$states, window = 100, step = 50,
                           max_lag = 20, min_stable = 3, threshold = 30)
  rc <- tdsReconfiguration(res, from = "A", to = "B")

  expect_s3_class(rc, "data.frame")
  expect_named(rc, c("from_node", "to_node", "tds_A", "tds_B", "delta"))
  expect_gt(rc$delta[1], 0)              # brain-heart strengthens A -> B
})

test_that("states with too few windows are dropped with a warning", {
  d <- make_two_state(n = 4000)
  st <- d$states
  st[1:5] <- "tiny"                      # a sliver of a third state
  expect_warning(
    tdsNetworkByState(d$X, st, window = 100, step = 50, max_lag = 20,
                      min_stable = 3, threshold = 30, min_windows = 3),
    "dropped")
})

test_that("states must align with the series length", {
  d <- make_two_state()
  expect_error(
    tdsNetworkByState(d$X, d$states[1:10], window = 100, threshold = 30),
    "one label per row")
})
