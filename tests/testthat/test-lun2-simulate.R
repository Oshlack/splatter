context("Lun2 simulation")

test_that("Lun2 simulation output is valid", {
    expect_true(validObject(lun2Simulate()))
    expect_true(validObject(lun2Simulate(de.nGenes = 100)))
    expect_true(validObject(lun2Simulate(zinb = TRUE)))
})

test_that("Gene sampling works", {
    expect_warning(lun2Simulate(nGenes = 10),
                   paste("Number of gene parameters does not equal nGenes.",
                         "Gene parameters will be sampled."))
})

test_that("Cell sampling works", {
    expect_warning(lun2Simulate(cell.plates = c(1, 1, 2)),
                   paste("Number of lib.sizes not equal to nCells.",
                         "lib.sizes will be sampled."))
})
