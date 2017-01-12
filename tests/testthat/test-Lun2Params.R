context("Lun2Params")

test_that("nCells checks work", {
    params <- newLun2Params()
    expect_error(setParam(params, "nCells", 1),
                 "nCells cannot be set directly, set cell.plates instead")
    expect_error(setParam(params, "nPlates", 1),
                 "nPlates cannot be set directly, set cell.plates instead")
})

test_that("cell sampling works", {
    params <- newLun2Params()
    expect_warning(setParam(params, "cell.plates", 1),
                   paste("nCells has been changed. cell.libSizes will be",
                         "sampled to length nCells"))
})

test_that("gene sampling works", {
    params <- newLun2Params()
    expect_warning(setParam(params, "nGenes", 1),
                   paste("nGenes has been changed. Gene parameter vectors will",
                         "be sampled to length new nGenes."))
})