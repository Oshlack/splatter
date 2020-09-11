context("Splat simulations")

test.params <- newSplatParams(nGenes = 100, batchCells = c(5, 5),
                              group.prob = c(0.5, 0.5), lib.scale = 0)


test_that("splatSimulate output is valid", {
    expect_true(validObject(splatSimulate(test.params, method = "single")))
    expect_true(validObject(splatSimulate(test.params, method = "groups")))
    expect_true(validObject(splatSimulate(test.params, method = "paths",
                                          path.from = c(0, 1))))
})

# Make mock vcf file
genotypes <- c(rep("0/0", 20), rep("0/1", 2), rep("1/1", 1))
nsnps <- 1e5
vcf <- data.frame(list(V1=rep(22, nsnps), V2=sample(1:5e7, nsnps)))
vcf[, c("V3","V4","V5","V6","V7","V8","V9")] <- NA
vcf[, c("V10","V11","V12","V13","V14","V15")] <- sample(
    genotypes, 6*nsnps, replace=TRUE)

# Set up for population simulation
pop.params <- newPopParams(eqtl.n = 0.5, nGenes = 100) 
pop <- popSimulate(vcf = vcf, popParams = pop.params)

params <- newSplatParams(batchCells = c(10, 10),
                         group.prob = c(0.5, 0.5), lib.scale = 0)

test_that("splatSimulatePop output is valid", {
    expect_true(validObject(splatSimulatePop(pop = pop$means, 
                                              params=params, 
                                              method = "single")))
    expect_true(validObject(splatSimulatePop(pop = pop$means, 
                                              params=params, 
                                              method = "groups")))
    expect_true(validObject(splatSimulatePop(pop = pop$means, 
                                              params=params, 
                                              method = "paths",
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


