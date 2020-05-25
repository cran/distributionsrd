context("Double-Pareto Lognormal functions")

test_that("Density, distribution function, quantile function, raw moments and random generation for the Double-Pareto Lognormal distribution work correctly", {
  expect_equal(2, qleftparetolognormal(pleftparetolognormal(2, log.p = TRUE), log.p = TRUE), tolerance = 1e-4)
  expect_equal(pleftparetolognormal(2), mleftparetolognormal(truncation = 2))

  x <- rleftparetolognormal(1e5)
  expect_equal(mean(x), mleftparetolognormal(r = 1, lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mleftparetolognormal(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE), tolerance = 1e-1)
})


test_that("Comparing probabilites of power-transformed variables", {
  coeff <- leftparetolognormal_plt(shape1 = 1.5, meanlog = -0.5, sdlog = 0.5, a = 5, b = 7)$coefficients
  expect_equal(
    pleftparetolognormal(3, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5),
    pleftparetolognormal(5 * 3^7, shape1 = coeff[["shape1"]], meanlog = coeff[["meanlog"]], sdlog = coeff[["sdlog"]])
  )

  coeff <- leftparetolognormal_plt(shape1 = 1.5, meanlog = -0.5, sdlog = 0.5, a = 5, b = 7, inv = TRUE)$coefficients
  expect_equal(
    pleftparetolognormal(5 * 0.9^7, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5),
    pleftparetolognormal(0.9, shape1 = coeff[["shape1"]], meanlog = coeff[["meanlog"]], sdlog = coeff[["sdlog"]])
  )
})
