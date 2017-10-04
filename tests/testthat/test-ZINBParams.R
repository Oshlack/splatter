context("ZINBParams")

test_that("constructor is valid", {
    expect_true(validObject(newZINBParams()))
})

test_that("nGenes checks work", {
    params <- newZINBParams()
    expect_error(setParam(params, "nGenes", 1),
                 "nGenes cannot be set directly")
    expect_error(setParam(params, "nCells", 1),
                 "nGenes cannot be set directly")
})
