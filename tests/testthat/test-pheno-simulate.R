context("PhenoPath simulation")

test_that("PhenoPath simulation output is valid", {
    expect_true(validObject(phenoSimulate()))
})
