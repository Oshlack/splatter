context("mfa simulation")

test_that("mfa simulation output is valid", {
    skip_if_not_installed("mfa")

    sim <- expect_warning(mfaSimulate(), "deprecated")
    expect_true(validObject(sim))

    sim <- expect_warning(mfaSimulate(dropout.present = TRUE), "deprecated")
    expect_true(validObject(sim))
})
