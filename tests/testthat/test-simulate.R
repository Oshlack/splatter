context("splatter simulations")

test_that("one group switches to single mode", {
  expect_warning(splat(method = "groups"),
                 "nGroups is 1, switching to single mode")
  expect_silent(splat(method = "paths", verbose = FALSE))
})
