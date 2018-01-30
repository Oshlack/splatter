context("PhenoParams")

params <- newPhenoParams()

test_that("constructor is valid", {
    expect_true(validObject(params))
})

test_that("printing works", {
    expect_output(show(params), "A Params object of class PhenoParams")
})

test_that("nGenes checks work", {
    expect_error(setParam(params, "nGenes", 1),
                 "nGenes cannot be set directly")
    pp <- setParam(params, "n.de", 0)
    total <- getParam(pp, "n.de") + getParam(pp, "n.pst") +
             getParam(pp, "n.pst.beta") + getParam(pp, "n.de.pst.beta")
    expect_equal(getParam(pp, "nGenes"), total)
})
