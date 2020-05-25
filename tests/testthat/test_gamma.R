context("Gamma functions")

test_that("Raw moments for the Gamma distribution work correctly", {
  expect_equal(pgamma(2, shape = 2, rate = 1), mgamma(truncation = 2))

  x <- rgamma(1e5, shape = 2, rate = 1)
  expect_equal(mean(x), mgamma(r = 1, lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mgamma(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE), tolerance = 1e-1)
})
