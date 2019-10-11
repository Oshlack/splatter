context("phenoEstimate")

library(scater)
set.seed(1)
counts <- counts(mockSCE())

test_that("phenoEstimate works", {
    skip_if_not_installed("phenopath")
    params <- phenoEstimate(counts)
    expect_true(validObject(params))
})
