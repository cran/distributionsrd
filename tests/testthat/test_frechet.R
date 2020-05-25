context("frechet functions")

test_that("Density, distribution function, quantile function, raw moments and random generation for the frechet distribution work correctly", {
  expect_equal(pfrechet(q = 5, shape = 2, scale = 1.5), 1 - pweibull(q = 1 / 5, shape = 2, scale = 1 / 1.5))

  expect_equal(2, qfrechet(pfrechet(2, log.p = T), log.p = T))
  expect_equal(pfrechet(2), mfrechet(truncation = 2))

  x <- rfrechet(1e5, scale = 1)
  expect_equal(mean(x), mfrechet(r = 1, lower.tail = FALSE, scale = 1), tolerance = 1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mfrechet(r = 1, truncation = quantile(x, 0.1), scale = 1, lower.tail = FALSE), tolerance = 1)
})


test_that("Comparing probabilites of power-transformed variables", {
  coeff <- frechet_plt(shape = 2, scale = 1, a = 5, b = 7)$coefficients
  expect_equal(pfrechet(3, shape = 2, scale = 1), pfrechet(5 * 3^7, shape = coeff[["shape"]], scale = coeff[["scale"]]))

  coeff <- frechet_plt(shape = 2, scale = 1, a = 5, b = 7, inv = TRUE)$coefficients
  expect_equal(pfrechet(5 * 0.8^7, shape = 2, scale = 1), pfrechet(0.8, shape = coeff[["shape"]], scale = coeff[["scale"]]))
})
