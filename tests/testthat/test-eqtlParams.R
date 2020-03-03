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
    expect_error(setParam(params, "bulkcv.param", data.frame(
        start = c(0, 10), shape = c(10, 2), rate = c(7, 3))),
                 "Need to set bulkcv.bins to length of bulkcv.param")

})


test_that("setParams order doesn't matter", {
    expect_silent(setParams(params, eqtl.n = 10, eqtl.maf = 0.5))
    expect_silent(setParams(params, eqtl.maf = 0.5, eqtl.n = 10))
})
