context("SplatParams")

test_that("nGroups checks work", {
  params <- newSplatParams()
  expect_error(setParam(params, "nCells", 1),
               "nCells cannot be set directly, set groupCells instead")
  expect_error(setParam(params, "nGroups", 1),
               "nGroups cannot be set directly, set groupCells instead")
})

test_that("path.from checks work", {
    params <- newSplatParams()
    params <- setParams(params, groupCells = c(10, 10))
    params <- setParamUnchecked(params, "path.from", c(0, 1))
    expect_silent(validObject(params))
    params <- setParamUnchecked(params, "path.from", c(0, 3))
    expect_error(validObject(params),
                 paste('invalid class "SplatParams" object: path.from:',
                       "All elements must be <= 2"))
    params <- setParamUnchecked(params, "path.from", c(1, 0))
    expect_error(validObject(params), "path cannot begin at itself")
    params <- newSplatParams()
    params <- setParams(params, groupCells = c(10, 10, 10))
    params <- setParamUnchecked(params, "path.from", c(2, 1, 1))
    expect_error(validObject(params), "origin must be specified in path.from")
})