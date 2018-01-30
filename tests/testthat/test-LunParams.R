context("LunParams")

params <- newLunParams()

test_that("printing works", {
    expect_output(show(params), "A Params object of class LunParams")
})

test_that("nCells checks work", {
    expect_error(setParam(params, "nCells", 1),
                 "nCells cannot be set directly, set groupCells instead")
})

test_that("nGroups checks work", {
    expect_error(setParam(params, "nGroups", 1),
                 "nGroups cannot be set directly, set groupCells instead")
})
