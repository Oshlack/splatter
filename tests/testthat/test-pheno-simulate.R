context("PhenoPath simulation")

test_that("PhenoPath simulation output is valid", {
    skip_if_not_installed("phenopath")
    expect_true(validObject(phenoSimulate()))
})
