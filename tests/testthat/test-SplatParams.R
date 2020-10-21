context("SplatParams")

params <- newSplatParams()

test_that("printing works", {
    expect_output(show(params), "A Params object of class SplatParams")
})

test_that("nBatches checks work", {
  expect_error(setParam(params, "nCells", 1),
               "nCells cannot be set directly, set batchCells instead")
  expect_error(setParam(params, "nBatches", 1),
               "nBatches cannot be set directly, set batchCells instead")
})

test_that("nGroups checks work", {
    expect_error(setParam(params, "nGroups", 1),
                 "nGroups cannot be set directly, set group.prob instead")
})


### These tests are also run in test-SplatPopParams.R, please update both
test_that("path.from checks work", {
    pp <- setParams(params, group.prob = c(0.5, 0.5))
    pp <- setParamUnchecked(pp, "path.from", c(0, 1))
    expect_silent(validObject(pp))
    pp <- setParamUnchecked(pp, "path.from", c(0, 3))
    expect_error(validObject(pp), "invalid class")
    pp <- setParamUnchecked(pp, "path.from", c(1, 0))
    expect_error(validObject(pp), "path cannot begin at itself")
    pp <- newSplatParams()
    pp <- setParams(pp, group.prob = c(0.3, 0.3, 0.4))
    pp <- setParamUnchecked(pp, "path.from", c(2, 1, 1))
    expect_error(validObject(pp), "origin must be specified in path.from")
    pp <- setParams(params, group.prob = c(0.5, 0.5), path.from = c(0, 1))
    expect_warning(setParam(pp, "group.prob", 1),
                   "nGroups has changed, resetting path.from")
    pp <- newSplatParams()
    pp <- setParams(pp, group.prob = c(0.25, 0.25, 0.25, 0.25))
    pp <- setParamUnchecked(pp, "path.from", c(0, 4, 2, 3))
    expect_error(validObject(pp), "path.from cannot contain cycles")
})

### These tests are also run in test-SplatPopParams.R, please update both
test_that("dropout.type checks work", {
    expect_error(setParam(params, "dropout.type", "cell"),
                 "dropout.type cannot be set to 'cell'")
    pp <- setParams(params, dropout.mid = rep(1, 100),
                    dropout.shape = rep(1, 100))
    expect_silent(setParam(pp, "dropout.type", "cell"))
    expect_error(setParam(params, "dropout.type", "a"),
                 "dropout.type must be one of: ")
})


test_that("setParams order doesn't matter", {
    expect_silent(setParams(params, group.prob = c(0.5, 0.5),
                            de.facLoc = c(0.1, 5)))
    expect_silent(setParams(params, de.facLoc = c(0.1, 5),
                            group.prob = c(0.5, 0.5)))
})
