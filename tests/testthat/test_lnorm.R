context("Lognormal functions")

test_that("Raw moments for the Lognormal distribution work correctly", {
  expect_equal(plnorm(2, meanlog = -0.5, sdlog = 0.5), mlnorm(truncation = 2))

  x <- rlnorm(1e5, meanlog = -0.5, sdlog = 0.5)
  expect_equal(mean(x), mlnorm(r = 1, lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mlnorm(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE), tolerance = 1e-1)
})


test_that("Comparing probabilites of power-transformed variables", {
  coeff <- lnorm_plt(meanlog = -0.5, sdlog = 0.5, a = 5, b = 7)$coefficients
  expect_equal(plnorm(3, meanlog = -0.5, sdlog = 0.5), plnorm(5 * 3^7, meanlog = coeff[["meanlog"]], sdlog = coeff[["sdlog"]]))

  coeff <- lnorm_plt(meanlog = -0.5, sdlog = 0.5, a = 5, b = 7, inv = TRUE)$coefficients
  expect_equal(plnorm(5 * 0.8^7, meanlog = -0.5, sdlog = 0.5), plnorm(0.8, meanlog = coeff[["meanlog"]], sdlog = coeff[["sdlog"]]))
})
