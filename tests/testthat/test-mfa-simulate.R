context("mfa simulation")

test_that("mfa simulation output is valid", {
    expect_true(validObject(mfaSimulate()))
    expect_true(validObject(mfaSimulate(dropout.present = TRUE)))
})
