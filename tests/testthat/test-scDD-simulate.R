context("scDD simulation")

test_that("scDD simulation output is valid", {
    params <- newSCDDParams(nDE = 5, nDP = 5, nDM = 5, nDB = 5, nEE = 5,
                            nEP = 5)
    expect_true(validObject(scDDSimulate(params)))
})