context("splatPopSimulate")

library(vcfR)
set.seed(1)

n.samples <- 6
n.genes <- 100
vcf <- mockVCF(n.samples=n.samples)

params <- setParams(newSplatPopParams(), eqtl.n = 10, nGenes=n.genes)

test_that("splatPopSimulate output is valid and works", {
    expect_true(validObject(splatPopSimulate(params = params, vcf = vcf)))
    pop <- splatPopSimulateMeans(params = params, vcf = vcf)
    expect_false(any(is.na(pop$means)))
    expect_false(any(sapply(pop$means, is.infinite)))
    expect_length(pop$means, n.samples)
    expect_length(pop$key, 13) # Number of columns expected in splatPop key
    expect_true(validObject(splatPopSimulateSC(sim.means = pop$means,
                                               key = pop$key,
                                               params = params)))
})


test_that("splatPopSimulate on multiple groups output is valid and works", {
    params.g2 <- setParams(params, group.prob = c(0.5, 0.5), nGenes=n.genes)
    eqtl <- splatPopSimulate(params = params.g2, vcf = vcf)
    expect_true(validObject(eqtl))

})

test_that("splatPopSimulate can read genes from gff data.frame", {
    gff <- mockGFF(n.genes = n.genes)
    sim <- splatPopSimulateMeans(vcf = vcf, params = params, gff=gff)
    expect_true(validObject(sim))
})
