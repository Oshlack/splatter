context("ZINB-WaVE simulation")

test_that("ZINB-WaVE simulation output is valid", {
    skip_if_not_installed("zinbwave")
    expect_true(validObject(zinbSimulate()))
})
