context("popSimulate")

set.seed(1)

# Random vcf matrix
genotypes <- c(rep("0/0", 20), rep("0/1", 2), rep("1/1", 1))
nsnps <- 1e5
vcf <- data.frame(list(V1=rep(22, nsnps), V2=sample(1:5e7, nsnps)))
vcf[, c("V3","V4","V5","V6","V7","V8","V9")] <- NA
vcf[, c("V10","V11","V12","V13","V14","V15")] <- sample(
    genotypes, 6*nsnps, replace=TRUE)

pop.params <- setParams(newPopParams(), eqtl.n = 10, seed = 42, nGenes=100)

test_that("popSimulate output is valid and works", {
    pop <- popSimulate(popParams = pop.params, vcf = vcf)
    expect_true(validObject(pop))
    expect_false(any(is.na(pop$means)))
    expect_false(any(sapply(pop$means, is.infinite)))
    expect_length(pop$means, 6)
})

pop.params.g2 <- setParams(newPopParams(), eqtl.n = 10, seed = 42, 
                           eqtl.groups = 2, nGenes=100)
test_that("popSimulate on multiple groups output is valid and works", {
    eqtl <- popSimulate(popParams = pop.params.g2, vcf = vcf)
    expect_true(validObject(eqtl))
    expect_false(any(is.na(eqtl$means[[1]])))
    expect_false(any(sapply(eqtl$means[[1]], is.infinite)))
    expect_length(eqtl$means, 2)
})


