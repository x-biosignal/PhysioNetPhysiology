test_that("physioNodeMatrix builds from a named list", {
  L <- list(a = 1:10, b = rnorm(10), c = runif(10))
  M <- physioNodeMatrix(L)
  expect_true(is.matrix(M))
  expect_equal(dim(M), c(10, 3))
  expect_equal(colnames(M), c("a", "b", "c"))
})

test_that("physioNodeMatrix accepts matrices and data frames", {
  df <- data.frame(a = 1:5, b = 6:10)
  M <- physioNodeMatrix(df)
  expect_equal(colnames(M), c("a", "b"))
  expect_equal(physioNodeMatrix(as.matrix(df))[, "b"], 6:10,
               ignore_attr = TRUE)
})

test_that("physioNodeMatrix rejects ragged or non-numeric input", {
  expect_error(physioNodeMatrix(list(a = 1:10, b = 1:5)), "same length")
  expect_error(physioNodeMatrix(list(a = 1:5, b = letters[1:5])), "numeric")
  expect_error(physioNodeMatrix(data.frame(a = 1:3, b = letters[1:3])),
               "numeric")
})

test_that("physioNodeMatrix can rename and standardize", {
  L <- list(a = rnorm(50, mean = 10), b = rnorm(50, mean = -3))
  M <- physioNodeMatrix(L, nodes = c("brain", "heart"), standardize = TRUE)
  expect_equal(colnames(M), c("brain", "heart"))
  expect_equal(colMeans(M), c(brain = 0, heart = 0), tolerance = 1e-8)
})
