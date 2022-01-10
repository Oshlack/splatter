context("lun2Estimate")

library(scater)
set.seed(1)
sce <- mockSCE(ngenes = 100)

test_that("lun2Estimate works", {
    plates <- as.numeric(factor(colData(sce)$Mutation_Status))
    params <- lun2Estimate(sce, plates, min.size = 20)
    expect_true(validObject(params))
})
