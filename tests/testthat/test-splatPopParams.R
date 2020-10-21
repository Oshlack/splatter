context("popParams")

params <- newSplatPopParams()

test_that("printing works", {
    expect_output(show(params), "A Params object of class SplatPopParams")
})

test_that("nCells checks work", {
    expect_error(setParam(params, "nCells", 1),
                   "nCells cannot be set directly, set batchCells instead")
})



test_that("CV params checks work", {
    expect_error(setParam(params, "pop.cv.param", data.frame(
        start = c(0, 10), shape = c(10, 2), rate = c(7, 3))),
                 "Need to set pop.cv.bins to length of pop.cv.param")

})

test_that("setParams order doesn't matter", {
    expect_silent(setParams(params, eqtl.n = 10, eqtl.maf.max = 0.4))
    expect_silent(setParams(params, eqtl.maf.max = 0.4, eqtl.n = 10))
})

### These tests are also run in test-splatParams.R, please update both
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

### These tests are also run in test-splatParams.R, please update both
test_that("dropout.type checks work", {
    expect_error(setParam(params, "dropout.type", "cell"),
                 "dropout.type cannot be set to 'cell'")
    pp <- setParams(params, dropout.mid = rep(1, 100),
                    dropout.shape = rep(1, 100))
    expect_silent(setParam(pp, "dropout.type", "cell"))
    expect_error(setParam(params, "dropout.type", "a"),
                 "dropout.type must be one of: ")
})
