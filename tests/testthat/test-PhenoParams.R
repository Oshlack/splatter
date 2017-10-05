context("PhenoParams")

test_that("constructor is valid", {
    expect_true(validObject(newPhenoParams()))
})

test_that("nGenes checks work", {
    params <- newPhenoParams()
    expect_error(setParam(params, "nGenes", 1),
                 "nGenes cannot be set directly")
    params <- setParam(params, "n.de", 0)
    total <- getParam(params, "n.de") + getParam(params, "n.pst") +
             getParam(params, "n.pst.beta") + getParam(params, "n.de.pst.beta")
    expect_equal(getParam(params, "nGenes"), total)
})
