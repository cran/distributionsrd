context("Composite functions")

test_that("Density, distribution function, quantile function, raw moments and random generation for the composite distribution work correctly", {

  ## Three-component distribution
  dist <- c("invpareto", "lnorm", "pareto")
  coeff <- c(coeff2.meanlog = -0.5, coeff2.sdlog = 0.5, coeff3.k = 1.5, coeff1.k = 1.5)


  # Demonstration of log functionality for probability and quantile function
  expect_equal(2, qcomposite(pcomposite(2, dist = dist, coeff = coeff, log.p = TRUE), dist = dist, coeff = coeff, log.p = TRUE))

  # The zeroth truncated moment is equivalent to the probability function
  expect_equal(pcomposite(2, dist = dist, coeff = coeff), mcomposite(truncation = 2, dist = dist, coeff = coeff))

  # The (truncated) first moment is equivalent to the mean of a (truncated) random sample, for large enough samples.
  coeff <- c(coeff2.meanlog = -0.5, coeff2.sdlog = 0.5, coeff3.k = 3, coeff1.k = 1.5)
  x <- rcomposite(1e5, dist = dist, coeff = coeff)

  expect_equal(mean(x), mcomposite(r = 1, lower.tail = FALSE, dist = dist, coeff = coeff), tolerance = 1e-1)

  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mcomposite(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE, dist = dist, coeff = coeff), tolerance = 1e-1)

  ## Two-component distribution
  dist <- c("lnorm", "pareto")
  coeff <- coeff <- c(coeff2.k = 1.5, coeff1.meanlog = -0.5, coeff1.sdlog = 0.5)

  # Demonstration of log functionality for probability and quantile function
  expect_equal(2, qcomposite(pcomposite(2, dist = dist, coeff = coeff, log.p = TRUE), dist = dist, coeff = coeff, log.p = TRUE))

  # The zeroth truncated moment is equivalent to the probability function
  expect_equal(pcomposite(2, dist = dist, coeff = coeff), mcomposite(truncation = 2, dist = dist, coeff = coeff))

  # The (truncated) first moment is equivalent to the mean of a (truncated) random sample, for large enough samples.
  coeff <- c(coeff1.meanlog = -0.5, coeff1.sdlog = 0.5, coeff2.k = 3)
  x <- rcomposite(1e5, dist = dist, coeff = coeff)

  expect_equal(mean(x), mcomposite(r = 1, lower.tail = FALSE, dist = dist, coeff = coeff), tolerance = 1e-1)

  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mcomposite(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE, dist = dist, coeff = coeff), tolerance = 1e-1)
})


test_that("Comparing probabilites of power-transformed variables", {

  ## Comparing probabilites of power-law transformed transformed variables
  dist <- c("invpareto", "lnorm", "pareto")
  coeff <- c(coeff2.meanlog = -0.5, coeff2.sdlog = 0.5, coeff3.k = 1.5, coeff1.k = 1.5)

  newcoeff <- composite_plt(dist = dist, coeff = coeff, a = 5, b = 7)$coefficients
  expect_equal(pcomposite(3, dist = dist, coeff = coeff), pcomposite(5 * 3^7, dist = dist, coeff = newcoeff))

  newcoeff <- composite_plt(dist = dist, coeff = coeff, a = 5, b = 3, inv = T)$coefficients
  expect_equal(pcomposite(5 * 0.9^3, dist = dist, coeff = coeff), pcomposite(0.9, dist = dist, coeff = newcoeff))
})
