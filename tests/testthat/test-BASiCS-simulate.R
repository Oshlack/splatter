context("BASiCS simulation")

test_that("BASiCS simulation output is valid", {
    skip_if_not_installed("BASiCS")
    expect_true(validObject(BASiCSSimulate()))
})

test_that("BASiCS simulation works when spike.means is sampled", {
    skip_if_not_installed("BASiCS")
    params <- newBASiCSParams(spike.means = c(12.93, 30.19, 1010.72, 7.9))
    expect_warning(sim <- BASiCSSimulate(params), "spike.means will be sampled")
    expect_true(validObject(sim))
})
