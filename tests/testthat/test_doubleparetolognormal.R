context("Double-Pareto Lognormal functions")

test_that("Density, distribution function, quantile function, raw moments and random generation for the Double-Pareto Lognormal distribution work correctly", {
  expect_equal(2, qdoubleparetolognormal(pdoubleparetolognormal(2, log.p = TRUE), log.p = TRUE), tolerance = 1e-4)
  expect_equal(pdoubleparetolognormal(2), mdoubleparetolognormal(truncation = 2))

  x <- rdoubleparetolognormal(1e5, shape2 = 3)
  expect_equal(mean(x), mdoubleparetolognormal(r = 1, shape2 = 3, lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mdoubleparetolognormal(r = 1, shape2 = 3, truncation = quantile(x, 0.1), lower.tail = FALSE), tolerance = 1e-1)
})


test_that("Comparing probabilites of power-transformed variables", {
  coeff <- doubleparetolognormal_plt(shape1 = 1.5, shape2 = 3, meanlog = -0.5, sdlog = 0.5, a = 5, b = 7)$coefficients
  expect_equal(
    pdoubleparetolognormal(3, shape1 = 1.5, shape2 = 3, meanlog = -0.5, sdlog = 0.5),
    pdoubleparetolognormal(5 * 3^7, shape1 = coeff[["shape1"]], shape2 = coeff[["shape2"]], meanlog = coeff[["meanlog"]], sdlog = coeff[["sdlog"]])
  )

  coeff <- doubleparetolognormal_plt(shape1 = 1.5, shape2 = 3, meanlog = -0.5, sdlog = 0.5, a = 5, b = 7, inv = TRUE)$coefficients
  expect_equal(
    pdoubleparetolognormal(5 * 0.9^7, shape1 = 1.5, shape2 = 3, meanlog = -0.5, sdlog = 0.5),
    pdoubleparetolognormal(0.9, shape1 = coeff[["shape1"]], shape2 = coeff[["shape2"]], meanlog = coeff[["meanlog"]], sdlog = coeff[["sdlog"]])
  )
})
