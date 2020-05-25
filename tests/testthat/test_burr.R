context("Burr functions")

test_that("Density, distribution function, quantile function, raw moments and random generation for the Burr distribution work correctly", {
  expect_equal(2, qburr(pburr(2, log.p = TRUE), log.p = TRUE))
  expect_equal(pburr(2), mburr(truncation = 2))

  x <- rburr(1e5, shape2 = 3)
  expect_equal(mean(x), mburr(r = 1, shape2 = 3, lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mburr(r = 1, shape2 = 3, truncation = as.numeric(quantile(x, 0.1)), lower.tail = FALSE), tolerance = 1e-1)
})

test_that("Comparing probabilites of power-transformed variables", {
  coeff <- burr_plt(shape1 = 2, shape2 = 3, scale = 1, a = 5, b = 7)$coefficients
  expect_equal(pburr(3, shape1 = 2, shape2 = 3, scale = 1), pburr(5 * 3^7, shape1 = coeff[["shape1"]], shape2 = coeff[["shape2"]], scale = coeff[["scale"]]))

  coeff <- burr_plt(shape1 = 2, shape2 = 3, scale = 1, a = 5, b = 7, inv = TRUE)$coefficients
  expect_equal(pburr(0.9, shape1 = coeff[["shape1"]], shape2 = coeff[["shape2"]], scale = coeff[["scale"]]), pburr(5 * 0.9^7, shape1 = 2, shape2 = 3, scale = 1))

  x <- rburr(1e5, shape1 = 2, shape2 = 3, scale = 1)
  coeff <- burr_plt(shape1 = 2, shape2 = 3, scale = 1, a = 2, b = 0.5)$coefficients
  y <- rburr(1e5, shape1 = coeff[["shape1"]], shape2 = coeff[["shape2"]], scale = coeff[["scale"]])

  expect_equal(mean(2 * x^0.5), mean(y), tolerance = 1e-1)
  expect_equal(mean(2 * x^0.5), mburr(r = 1, shape1 = coeff[["shape1"]], shape2 = coeff[["shape2"]], scale = coeff[["scale"]], lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(mean(y), mburr(r = 1, shape1 = coeff[["shape1"]], shape2 = coeff[["shape2"]], scale = coeff[["scale"]], lower.tail = FALSE), tolerance = 1e-1)
})
