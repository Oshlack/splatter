context("mfa simulation")

test_that("mfa simulation output is valid", {
    skip_if_not_installed("mfa")
    expect_true(validObject(mfaSimulate()))
    expect_true(validObject(mfaSimulate(dropout.present = TRUE)))
})
