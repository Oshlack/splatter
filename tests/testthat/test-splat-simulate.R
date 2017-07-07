context("Splatter simulations")

test.params <- newSplatParams(nGenes = 100, batchCells = c(5, 5),
                              group.prob = c(0.5, 0.5), lib.scale = 0)

test_that("splatSimulate output is valid", {
    expect_true(validObject(splatSimulate(test.params, method = "single")))
    expect_true(validObject(splatSimulate(test.params, method = "groups")))
    expect_true(validObject(splatSimulate(test.params, method = "paths")))
})

test_that("one group switches to single mode", {
    expect_warning(splatSimulate(test.params, method = "groups",
                                 group.prob = c(1)),
                   "nGroups is 1, switching to single mode")
    expect_silent(splatSimulate(test.params, method = "paths",
                                group.prob = c(1), verbose = FALSE))
})

test_that("infinite bcv.df is detected", {
    expect_warning(splatSimulate(test.params, bcv.df = Inf),
                   "'bcv.df' is infinite. This parameter will be ignored.")
})
