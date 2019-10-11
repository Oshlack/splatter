context("Lun simulation")

test_that("Lun simulation output is valid", {
    expect_true(validObject(lunSimulate()))
    expect_true(validObject(lunSimulate(groupCells = c(10, 10))))
})
