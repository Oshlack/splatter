context("SCDDParams")

test_that("constructor is valid", {
    expect_true(validObject(newSCDDParams()))
})

test_that("nGenes checks work", {
    params <- newSCDDParams()
    expect_error(setParam(params, "nGenes", 1),
                 paste("nGenes cannot be set directly, set nDE, nDP, nDM, nDB,",
                       "nEE or nEP instead"))
    params <- setParam(params, "nEE", 0)
    total <- getParam(params, "nDE") + getParam(params, "nDP") +
             getParam(params, "nDM") + getParam(params, "nDP") +
             getParam(params, "nEE") + getParam(params, "nEP")
    expect_equal(getParam(params, "nGenes"), total)
})