context("splatPopSimulate")

if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
    requireNamespace("preprocessCore", quietly = TRUE)) {
    set.seed(1)

    n.samples <- 6
    n.genes <- 10
    vcf <- mockVCF(n.samples = n.samples, n.snps = 1000)
    gff <- mockGFF(n.genes = n.genes)

    params <- newSplatPopParams(eqtl.n = 10, nGenes = n.genes, seed = 1)
}

test_that("splatPopSimulate output is valid and works", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    expect_true(validObject(splatPopSimulate(params = params, vcf = vcf)))
    pop <- splatPopSimulateMeans(params = params, vcf = vcf)
    expect_false(any(is.na(pop$means)))
    expect_false(any(sapply(pop$means, is.infinite)))
    expect_equal(ncol(pop$means), n.samples)
    expect_true(validObject(splatPopSimulateSC(
        sim.means = pop$means,
        key = pop$key,
        params = params
    )))
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
    sim <- splatPopSimulateMeans(vcf = vcf, params = params, gff = gff)
    expect_true(validObject(sim))
})

test_that("splatPopSimulate can simulate from empirical data directly", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    emp <- mockEmpiricalSet(seed = 1)

    sim <- splatPopSimulateMeans(
        vcf = emp$vcf, gff = emp$gff, eqtl = emp$eqtl,
        means = emp$means
    )

    expect_true(validObject(sim))
})

test_that("splatPopSimulate with nCells.sample gives different cell counts", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")

    params.nCells <- newSplatPopParams(
        nCells.sample = TRUE, nCells.shape = 1.5,
        nCells.rate = 0.01
    )

    sim <- splatPopSimulate(vcf = vcf, gff = gff, params = params.nCells)
    cell.counts <- table(colData(sim)$Sample)

    expect_false(all(cell.counts == cell.counts[1]))
})

test_that("splatPopSimulate seeds are reproducible", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")

    sim1 <- splatPopSimulate(
        vcf = vcf, gff = gff, params = params,
        batchCells = c(50, 50)
    )
    sim2 <- splatPopSimulate(
        vcf = vcf, gff = gff, params = params,
        batchCells = c(50, 50)
    )

    expect_identical(sim1, sim2)
})
