context("Double-Pareto Lognormal functions")

test_that("Density, distribution function, quantile function, raw moments and random generation for the Double-Pareto Lognormal distribution work correctly", {
  expect_equal(2, qrightparetolognormal(prightparetolognormal(2, log.p = TRUE), log.p = TRUE), tolerance = 1e-4)
  expect_equal(prightparetolognormal(2), mrightparetolognormal(truncation = 2))

  x <- rrightparetolognormal(1e5, shape2 = 3)
  expect_equal(mean(x), mrightparetolognormal(r = 1, shape2 = 3, lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mrightparetolognormal(r = 1, shape2 = 3, truncation = quantile(x, 0.1), lower.tail = FALSE), tolerance = 1e-1)
})


test_that("Comparing probabilites of power-transformed variables", {
  coeff <- rightparetolognormal_plt(shape2 = 3, meanlog = -0.5, sdlog = 0.5, a = 5, b = 7)$coefficients
  expect_equal(
    prightparetolognormal(3, shape2 = 3, meanlog = -0.5, sdlog = 0.5),
    prightparetolognormal(5 * 3^7, shape2 = coeff[["shape2"]], meanlog = coeff[["meanlog"]], sdlog = coeff[["sdlog"]])
  )

  coeff <- rightparetolognormal_plt(shape2 = 3, meanlog = -0.5, sdlog = 0.5, a = 5, b = 7, inv = TRUE)$coefficients
  expect_equal(
    prightparetolognormal(5 * 0.9^7, shape2 = 3, meanlog = -0.5, sdlog = 0.5),
    prightparetolognormal(0.9, shape2 = coeff[["shape2"]], meanlog = coeff[["meanlog"]], sdlog = coeff[["sdlog"]])
  )
})
