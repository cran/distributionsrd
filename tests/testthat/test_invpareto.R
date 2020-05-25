context("invpareto functions")

test_that("Density, distribution function, quantile function, raw moments and random generation for the invpareto distribution work correctly", {
  expect_equal(2, qinvpareto(pinvpareto(2, log.p = T), log.p = T))
  expect_equal(pinvpareto(2), minvpareto(truncation = 2))

  x <- rinvpareto(1e6)
  expect_equal(mean(x), minvpareto(r = 1, lower.tail = F), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), minvpareto(r = 1, truncation = quantile(x, 0.1), lower.tail = F), tolerance = 1e-1)
})


test_that("Comparing probabilites of power-transformed variables", {
  coeff <- invpareto_plt(xmax = 5, k = 2, a = 5, b = 7)$coefficients
  expect_equal(pinvpareto(3, k = 2, xmax = 5), pinvpareto(5 * 3^7, k = coeff[["k"]], xmax = coeff[["xmax"]]))

  coeff <- invpareto_plt(xmax = 2, k = 5, a = 5, b = 7, inv = TRUE)$coefficients
  expect_equal(pinvpareto(5 * 0.9^7, k = 5, xmax = 2), pinvpareto(0.9, k = coeff[["k"]], xmax = coeff[["xmax"]]))
})
