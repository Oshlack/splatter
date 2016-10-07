library(splatter)
context("splatParams object")

test_that("checkParams checks class", {
    expect_error(checkParams("a"),
                 "params does not belong to the splatParams class")
    expect_error(checkParams(1),
                 "params does not belong to the splatParams class")
    expect_error(checkParams(list()),
                 "params does not belong to the splatParams class")
})

test_that("checkParams checks numeric", {
    params <- splatParams()
    params$nGenes <- "A"
    expect_error(checkParams(params), "nGenes must be numeric")
    params$nGenes <- TRUE
    expect_error(checkParams(params), "nGenes must be numeric")
    params$nGenes <- 100
    expect_silent(checkParams(params))
})

test_that("checkParams checks integer", {
    params <- splatParams()
    params$nGenes <- 1.5
    expect_error(checkParams(params), "nGenes must be an integer")
    params$nGenes <- NA
    params$path$length <- 1.5
    expect_error(checkParams(params), "path.length must be an integer")
    params$path$length <- NA
    params$nGenes <- 100
    expect_silent(checkParams(params))
})

test_that("checkParams checks positive", {
    params <- splatParams()
    params$nGenes <- -1
    expect_error(checkParams(params), "nGenes must be positive")
    params$nGenes <- NA
    params$mean$rate <- -1
    expect_error(checkParams(params), "mean.rate must be positive")
    params$mean$rate <- 1
    expect_silent(checkParams(params))
})

test_that("checkParams checks prob", {
    params <- splatParams()
    params$out$prob <- 1.2
    expect_error(checkParams(params), "out.prob must be in the range 0-1")
    params$out$prob <- 0.5
    expect_silent(checkParams(params))
})

test_that("checkParams checks logical", {
    params <- splatParams()
    params$dropout$present <- "A"
    #expect_error(checkParams(params),
    #             "dropout.present must be logical (TRUE/FALSE)")
    params$dropout$present <- TRUE
    expect_silent(checkParams(params))
    params$dropout$present <- FALSE
    expect_silent(checkParams(params))
})

test_that("checkParams checks vectors allowed", {
    params <- splatParams()
    params$nGenes <- c(1, 2)
    expect_error(checkParams(params), "nGenes should be a single value")
    params$nGenes <- NA
    params$groupCells <- c(100, 200)
})

test_that("checkParams checks vector length", {
    params <- splatParams()
    params <- setParams(params, groupCells = c(100, 200))
    params$path$length <- 100
    expect_silent(checkParams(params))
    params$path$length <- c(100, 200, 300)
    expect_error(checkParams(params),
                 paste("length of path.length must be 1 or the length of the",
                       "groupCells parameter"))
    params$path$length <- c(100, 200)
    expect_silent(checkParams(params))
    params$groupCells <- NA
    expect_error(checkParams(params),
                 paste("length of path.length must be 1 or the length of the",
                       "groupCells parameter"))
})

test_that("checkParams checks vector is not NA", {
    params <- splatParams()
    params <- setParams(params, groupCells = c(100, 200))
    expect_silent(checkParams(params))
    params$groupCells <- c(100, NA)
    expect_error(checkParams(params),
                 "groupCells is a vector and contains NA values")
})

test_that("checkParams checks groupCells", {
    params <- defaultParams()
    params$nCells <- 1
    expect_error(checkParams(params),
                 "nCells, nGroups and groupCells are not consistent")
    params <- defaultParams()
    params$nGroups <- 10
    expect_error(checkParams(params),
                 "nCells, nGroups and groupCells are not consistent")
})

test_that("setParams sets correctly", {
    params <- splatParams()
    params <- setParams(params, nGenes = 100)
    expect_equal(params$nGenes, 100)
    params <- setParams(params, mean.rate = 0.5)
    expect_equal(params$mean$rate, 0.5)
    params <- setParams(params, groupCells = c(100, 200))
    expect_equal(params$groupCells, c(100, 200))
    params <- setParams(params, dropout.present = TRUE)
    expect_equal(params$dropout$present, TRUE)
})

test_that("setParmas sets groupCells correctly", {
    params <- splatParams()
    expect_error(setParams(params, nCells = 100),
                 "nCells cannot be set directly, set groupCells instead")
    expect_error(setParams(params, nGroups = 10),
                 "nGroups cannot be set directly, set groupCells instead")
    params <- setParams(params, groupCells = 100)
    params <- setParams(params, groupCells = c(100, 200, 300))
})

test_that("getParams gets correctly", {
    params <- defaultParams()
    expect_equal(getParams(params, "nGenes"), c(nGenes = 10000))
    expect_equal(getParams(params, "mean.rate"), c(mean.rate = 0.3))
    expect_equal(getParams(params, c("nGenes", "mean.rate")),
                 c(nGenes = 10000, mean.rate = 0.3))
    params <- setParams(params, groupCells = c(100, 200))
    expect_equal(getParams(params, "groupCells"), c(100, 200))
    expect_equal(getParams(params, c("nGenes", "mean.rate", "groupCells")),
                 list(nGenes = 10000, mean.rate = 0.3,
                      groupCells = c(100, 200)))
})

test_that("mergeParams merges correctly", {
    params1 <- splatParams(nGenes = 100, mean.rate = 0.5,
                           groupCells = c(100, 200))
    params2 <- defaultParams()
    params <- mergeParams(params1, params2)
    expect_equal(params$nGenes, params1$nGenes)
    expect_equal(params$mean$rate, params1$mean$rate)
    expect_equal(params$groupCells, params1$groupCells)
    expect_equal(params$out$prob, params2$out$prob)
    expect_equal(params$dropout$present, params2$dropout$present)
})

test_that("constructor is valid", {
    expect_silent(checkParams(splatParams()))
})

test_that("defaultParams is valid", {
    expect_silent(checkParams(defaultParams()))
})