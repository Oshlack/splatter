context("Lun simulation")

test_that("Lun simulation output is valid", {
    expect_true(validObject(simSimple()))
    expect_true(validObject(simSimple(groupCells = c(10, 10))))
})