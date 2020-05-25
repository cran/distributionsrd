context("Weibull functions")

test_that("Raw moments for the Weibull distribution work correctly", {
  expect_equal(pweibull(2, shape = 2, scale = 1), mweibull(truncation = 2))

  x <- rweibull(1e5, shape = 2, scale = 1)
  expect_equal(mean(x), mweibull(r = 1, lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mweibull(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE), tolerance = 1e-1)
})


test_that("Comparing probabilites of power-transformed variables", {
  coeff <- weibull_plt(shape = 2, scale = 1, a = 5, b = 7)$coefficients
  expect_equal(pweibull(3, shape = 2, scale = 1), pweibull(5 * 3^7, shape = coeff[["shape"]], scale = coeff[["scale"]]))

  coeff <- weibull_plt(shape = 2, scale = 1, a = 5, b = 7, inv = TRUE)$coefficients
  expect_equal(pweibull(5 * 0.8^7, shape = 2, scale = 1), pweibull(0.8, shape = coeff[["shape"]], scale = coeff[["scale"]]))
})
