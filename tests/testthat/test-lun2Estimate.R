context("lun2Estimate")

library(scuttle)
set.seed(1)
sce <- mockSCE(ngenes = 100)

test_that("lun2Estimate works", {
    plates <- as.numeric(factor(colData(sce)$Mutation_Status))
    params <- expect_warning(
        lun2Estimate(sce, plates, min.size = 20),
        "more singular values"
    )
    expect_true(validObject(params))
})
