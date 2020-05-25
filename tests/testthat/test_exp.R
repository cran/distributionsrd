context("Exponential functions")

test_that("Raw moments for the exponential distribution work correctly", {
  expect_equal(pexp(2, rate = 1), mexp(truncation = 2))

  x <- rexp(1e5, rate = 1)
  expect_equal(mean(x), mexp(r = 1, lower.tail = FALSE), tolerance = 1e-1)
  expect_equal(sum(x[x > quantile(x, 0.1)]) / length(x), mexp(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE), tolerance = 1e-1)
})
