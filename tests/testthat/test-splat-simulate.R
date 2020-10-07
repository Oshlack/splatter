context("Splat simulations")

test.params <- newSplatParams(nGenes = 100, batchCells = c(5, 5),
                              group.prob = c(0.5, 0.5), lib.scale = 0)

test_that("splatSimulate output is valid", {
    expect_true(validObject(splatSimulate(test.params, method = "single")))
    expect_true(validObject(splatSimulate(test.params, method = "groups")))
    expect_true(validObject(splatSimulate(test.params, method = "paths",
                                          path.from = c(0, 1))))
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

test_that("dropout.type checks work", {
    pp <- setParams(test.params, dropout.type = "experiment")
    expect_true(validObject(splatSimulate(pp, method = "single")))
    pp <- setParams(pp, dropout.mid = 1:2)
    expect_error(splatSimulate(pp), "aren't length 1")
    pp <- setParams(test.params, group.prob = c(0.5, 0.5),
                    dropout.mid = c(1, 2), dropout.shape = c(-1, -0.5),
                    dropout.type = "group")
    expect_error(splatSimulate(pp), "groups have not been simulated")
})
