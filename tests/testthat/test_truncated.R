context("Truncated functions")

test_that("Density, distribution function, quantile function, raw moments and random generation for truncated distributions work correctly", {

  # Compare quantities
  expect_equal(dtruncdist(0.5), dlnorm(0.5))

  expect_equal(ptruncdist(2), plnorm(2))

  expect_equal(qtruncdist(0.25), qlnorm(0.25))

  expect_equal(mtruncdist(r = 0, truncation = 2), mlnorm(r = 0, truncation = 2, meanlog = 0, sdlog = 1))

  expect_equal(mtruncdist(r = 1, truncation = 2), mlnorm(r = 1, truncation = 2, meanlog = 0, sdlog = 1))
})
