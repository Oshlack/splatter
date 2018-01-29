context("simple simulation")

test_that("simple simulation output is valid", {
  expect_true(validObject(simpleSimulate()))
})
