context("SplatParams")

test_that("nBatches checks work", {
  params <- newSplatParams()
  expect_error(setParam(params, "nCells", 1),
               "nCells cannot be set directly, set batchCells instead")
  expect_error(setParam(params, "nBatches", 1),
               "nBatches cannot be set directly, set batchCells instead")
})

test_that("nGroups checks work", {
    params <- newSplatParams()
    expect_error(setParam(params, "nGroups", 1),
                 "nGroups cannot be set directly, set group.prob instead")
})

test_that("path.from checks work", {
    params <- newSplatParams()
    params <- setParams(params, group.prob = c(0.5, 0.5))
    params <- setParamUnchecked(params, "path.from", c(0, 1))
    expect_silent(validObject(params))
    params <- setParamUnchecked(params, "path.from", c(0, 3))
    expect_error(validObject(params),
                 paste('invalid class “SplatParams” object: path.from:',
                       "All elements must be <= 2"))
    params <- setParamUnchecked(params, "path.from", c(1, 0))
    expect_error(validObject(params), "path cannot begin at itself")
    params <- newSplatParams()
    params <- setParams(params, group.prob = c(0.3, 0.3, 0.4))
    params <- setParamUnchecked(params, "path.from", c(2, 1, 1))
    expect_error(validObject(params), "origin must be specified in path.from")
})
