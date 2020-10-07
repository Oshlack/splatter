context("splatPopSimulate")

library(vcfR)
set.seed(42)

n_samples <- 6
n_genes <- 100
vcf <- mock_vcf(n_samples=n_samples)

params <- setParams(newSplatPopParams(), eqtl.n = 10, seed = 42, nGenes=n_genes)

test_that("splatPopSimulate output is valid and works", {
    pop <- splatPopSimulate(params = params, vcf = vcf)
    expect_true(validObject(pop))
    expect_false(any(is.na(pop$means)))
    expect_false(any(sapply(pop$means, is.infinite)))
    expect_length(pop$means, n_samples)
    expect_length(pop$key, 12) # Number of columns expected in splatPop key
})

params.g2 <- setParams(params, group.prob = c(0.5, 0.5), nGenes=n_genes)
test_that("splatPopSimulate on multiple groups output is valid and works", {
    eqtl <- splatPopSimulate(params = params.g2, vcf = vcf)
    expect_true(validObject(eqtl))
    expect_false(any(is.na(eqtl$means[[1]])))
    expect_false(any(sapply(eqtl$means[[1]], is.infinite)))
    expect_length(eqtl$means, 2)
})

test_that("splatPopSimulate can read genes from gff data.frame", {
    gff <- mock_gff(n_genes = n_genes)
    params_gff <- setParams(params, random.genes = FALSE)
    sim <- splatPopSimulateMeans(vcf = vcf, params = params_gff, gff=gff)
    expect_true(validObject(sim))
})
