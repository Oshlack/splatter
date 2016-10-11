context("splatter simulations")

test.params <- defaultParams()
test.params <- setParams(test.params, nGenes = 100, groupCells = c(5, 5),
                         lib.scale = 0)

test_that("splat output is valid", {
    expect_true(validObject(splat(test.params, method = "single")))
    expect_true(validObject(splat(test.params, method = "groups")))
    expect_true(validObject(splat(test.params, method = "paths")))
})

test_that("one group switches to single mode", {
  expect_warning(splat(test.params, method = "groups", groupCells = c(10)),
                 "nGroups is 1, switching to single mode")
  expect_silent(splat(test.params, method = "paths", groupCells = c(10),
                      verbose = FALSE))
})
