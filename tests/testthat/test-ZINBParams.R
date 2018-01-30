context("ZINBParams")

params <- newZINBParams()

test_that("constructor is valid", {
    expect_true(validObject(params))
})

test_that("printing works", {
    expect_output(show(params), "A Params object of class ZINBParams")
})

test_that("nGenes checks work", {
    expect_error(setParam(params, "nGenes", 1),
                 "nGenes cannot be set directly")
    expect_error(setParam(params, "nCells", 1),
                 "nCells cannot be set directly")
})
