context("SCDDParams")

params <- newSCDDParams()

test_that("constructor is valid", {
    expect_true(validObject(params))
})

test_that("printing works", {
    expect_output(show(params), "A Params object of class SCDDParams")
})

test_that("nGenes checks work", {
    expect_error(setParam(params, "nGenes", 1),
                 paste("nGenes cannot be set directly, set nDE, nDP, nDM, nDB,",
                       "nEE or nEP instead"))
    pp <- setParam(params, "nEE", 0)
    total <- getParam(pp, "nDE") + getParam(pp, "nDP") +
             getParam(pp, "nDM") + getParam(pp, "nDP") +
             getParam(pp, "nEE") + getParam(pp, "nEP")
    expect_equal(getParam(pp, "nGenes"), total)
})
