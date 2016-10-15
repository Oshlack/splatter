context("Lun2 simulation")

test_that("Lun2 simulation output is valid", {
    expect_true(validObject(lun2Simulate()))
    expect_true(validObject(lun2Simulate(de.nGenes = 100)))
    expect_true(validObject(lun2Simulate(zinb = TRUE)))
})