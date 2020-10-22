context("splatPopSimulate")

if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
    requireNamespace("preprocessCore", quietly = TRUE)) {

    set.seed(1)

    n.samples <- 6
    n.genes <- 10
    vcf <- mockVCF(n.samples = n.samples, n.snps = 1000)

    params <- setParams(newSplatPopParams(), eqtl.n = 10, nGenes = n.genes)
}

test_that("splatPopSimulate output is valid and works", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    expect_true(validObject(splatPopSimulate(params = params, vcf = vcf)))
    pop <- splatPopSimulateMeans(params = params, vcf = vcf)
    expect_false(any(is.na(pop$means)))
    expect_false(any(sapply(pop$means, is.infinite)))
    expect_equal(ncol(pop$means), n.samples)
    expect_length(pop$key, 15) # Number of columns expected in splatPop key
    expect_true(validObject(splatPopSimulateSC(sim.means = pop$means,
                                               key = pop$key,
                                               params = params)))
})


test_that("splatPopSimulate on multiple groups output is valid and works", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    params.g2 <- setParams(params, group.prob = c(0.5, 0.5), nGenes = n.genes)
    eqtl <- splatPopSimulate(params = params.g2, vcf = vcf)
    expect_true(validObject(eqtl))

})

test_that("splatPopSimulate can read genes from gff data.frame", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    gff <- mockGFF(n.genes = n.genes)
    sim <- splatPopSimulateMeans(vcf = vcf, params = params, gff = gff)
    expect_true(validObject(sim))
})
