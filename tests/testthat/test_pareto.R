context("Pareto functions")

test_that("Density, distribution function, quantile function, raw moments and random generation for the Pareto distribution work correctly", {
  expect_equal(ppareto(q = 5, k = 2, xmin = 3), pexp(q = log(5 / 3), rate = 2))
  expect_equal(2, qpareto(ppareto(2, log.p = T), log.p = T))
  expect_equal(ppareto(2), mpareto(truncation = 2))

  x <- rpareto(1e6)
  expect_equal(mean(x), mpareto(r = 1, lower.tail = F), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mpareto(r = 1, truncation = quantile(x, 0.1), lower.tail = F), tolerance = 1e-1)
})


test_that("Comparing probabilites of power-transformed variables", {
  coeff <- pareto_plt(xmin = 2, k = 2, a = 5, b = 7)$coefficients
  expect_equal(ppareto(3, k = 2, xmin = 2), ppareto(5 * 3^7, k = coeff[["k"]], xmin = coeff[["xmin"]]))

  coeff <- pareto_plt(xmin = 2, k = 2, a = 5, b = 7, inv = TRUE)$coefficients
  expect_equal(ppareto(5 * 0.9^7, k = 2, xmin = 2), ppareto(0.9, k = coeff[["k"]], xmin = coeff[["xmin"]]))

  x <- rpareto(1e5, k = 2, xmin = 2)
  coeff <- pareto_plt(xmin = 2, k = 2, a = 2, b = 0.5)$coefficients
  y <- rpareto(1e5, k = coeff[["k"]], xmin = coeff[["xmin"]])

  expect_equal(mean(2 * x^0.5), mean(y), tolerance = 1e-1)
  expect_equal(mean(2 * x^0.5), mpareto(r = 1, k = coeff[["k"]], xmin = coeff[["xmin"]], lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(mean(y), mpareto(r = 1, k = coeff[["k"]], xmin = coeff[["xmin"]], lower.tail = FALSE), tolerance = 1e-1)
})
