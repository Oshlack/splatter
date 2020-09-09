context("eQTLParams")

params <- neweQTLParams()

test_that("printing works", {
    expect_output(show(params), "A Params object of class eQTLParams")
})


test_that("nCells checks work", {
    expect_warning(setParam(params, "nCells", 1),
                 "nCells parameter does not impact eQTL simulation, set nCells
             in the sc simulation Params object instead.")
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
