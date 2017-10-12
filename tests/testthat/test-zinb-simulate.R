context("ZINB-WaVE simulation")

test_that("ZINB-WaVE simulation output is valid", {
    expect_true(validObject(zinbSimulate()))
})
