context("Empirical functions")

test_that("Density, distribution function, quantile function, raw moments for the empirical distribution work correctly", {
  x <- rlnorm(1e5, meanlog = -0.5, sdlog = 0.5)


  # Compare empirical and parametric quantities
  expect_equal(dlnorm(0.5, meanlog = -0.5, sdlog = 0.5), dempirical(0.5, data = x), tolerance = 1e-1)

  expect_equal(plnorm(0.5, meanlog = -0.5, sdlog = 0.5), pempirical(0.5, data = x), tolerance = 1e-1)

  expect_equal(qlnorm(0.5, meanlog = -0.5, sdlog = 0.5), qempirical(0.5, data = x), tolerance = 1e-1)

  expect_equal(mlnorm(r = 0, truncation = 0.5, meanlog = -0.5, sdlog = 0.5), mempirical(r = 0, truncation = 0.5, data = x), tolerance = 1e-1)

  expect_equal(mlnorm(r = 1, truncation = 0.5, meanlog = -0.5, sdlog = 0.5), mempirical(r = 1, truncation = 0.5, data = x), tolerance = 1e-1)

  ## Demonstration of log functionailty for probability and quantile function
  expect_equal(as.numeric(quantile(x, 0.5, type = 1)), qempirical(p = pempirical(q = quantile(x, 0.5, type = 1), data = x, log.p = TRUE), data = x, log.p = TRUE), tolerance = 1e-1)

  ## The zeroth truncated moment is equivalent to the probability function
  expect_equal(pempirical(q = quantile(x, 0.5, type = 1), data = x), mempirical(truncation = quantile(x, 0.5, type = 1), data = x), tolerance = 1e-1)

  ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample, for large enough samples.
  expect_equal(mean(x), mempirical(r = 1, data = x, truncation = 0, lower.tail = FALSE), tolerance = 1e-1)

  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mempirical(r = 1, data = x, truncation = quantile(x, 0.1), lower.tail = FALSE), tolerance = 1e-1)
})
