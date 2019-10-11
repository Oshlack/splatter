context("simpleEstimate")

library(scater)
set.seed(1)
counts <- counts(mockSCE())

test_that("simpleEstimate works", {
    params <- simpleEstimate(counts)
    expect_true(validObject(params))
})

test_that("simpleEstimate works with SingleCellExperiment", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts)
    )
    params <- simpleEstimate(sce)
    expect_true(validObject(params))
})

test_that("simpleEstimate works with SingleCellExperiment without counts", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(TEST = counts)
    )
    expect_warning(simpleEstimate(sce), "counts assay is missing")
})
