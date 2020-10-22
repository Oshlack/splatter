context("splatPopParams")

if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
    requireNamespace("preprocessCore", quietly = TRUE)) {
    params <- newSplatPopParams()
}

test_that("printing works", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    expect_output(show(params), "A Params object of class SplatPopParams")
})

test_that("nCells checks work", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    expect_error(setParam(params, "nCells", 1),
                   "nCells cannot be set directly, set batchCells instead")
})

test_that("CV params checks work", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    expect_error(setParam(params, "pop.cv.param", data.frame(
        start = c(0, 10), shape = c(10, 2), rate = c(7, 3))),
                 "Need to set pop.cv.bins to length of pop.cv.param")

})

test_that("setParams order doesn't matter", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    expect_silent(setParams(params, eqtl.n = 10, eqtl.maf.max = 0.4))
    expect_silent(setParams(params, eqtl.maf.max = 0.4, eqtl.n = 10))
})
