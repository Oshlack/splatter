context("phenoEstimate")

library(scater)
data("sc_example_counts")

test_that("phenoEstimate works", {
    skip_if_not_installed("phenopath")
    params <- phenoEstimate(sc_example_counts)
    expect_true(validObject(params))
})
