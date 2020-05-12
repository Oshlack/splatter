context("eQTLSimulate")

data(ex_snps)
ex_snps <- ex_snps[, 1:19]
eparams <- setParams(neweQTLParams(), list(eqtl.n = 10, seed = 42))

test_that("eQTLSimulate output is valid and works", {
    eqtl <- eQTLSimulate(eQTLparams = eparams, vcf = ex_snps)
    expect_true(validObject(eqtl))
    expect_false(any(is.na(eqtl$means)))
    expect_false(any(sapply(eqtl$means, is.infinite)))
    expect_length(eqtl$means, 10)
})

eparams.g2 <- setParams(neweQTLParams(), list(eqtl.n = 10, seed = 42, 
                                              eqtl.groups = 2))
test_that("eQTLSimulate on multiple groups output is valid and works", {
    eqtl <- eQTLSimulate(eQTLparams = eparams.g2, vcf = ex_snps)
    expect_true(validObject(eqtl))
    expect_false(any(is.na(eqtl$means[[1]])))
    expect_false(any(sapply(eqtl$means[[1]], is.infinite)))
    expect_length(eqtl$means, 2)
})

test_that("input data checks", {
    gff_bad <- ex_gff
    gff_bad$V3 <- 'na'
    expect_error(eQTLSimulate(gff = gff_bad),
                 "GFF file did not contain gene features or other issue with 
            file format. See example data.")
    
    vcf_bad <- ex_snps
    vcf_bad$V1 <- NULL
    expect_error(eQTLSimulate(vcf = vcf_bad, eQTLparams = eparams),
                 "snps not in the expected VCF format. See example data.")
    
    vcf_bad2 <- ex_snps[1:5,]
    expect_error(eQTLSimulate(vcf = vcf_bad2), 
    "Not enough SNPs within desired MAF range. Increase the
                    eqtl.mafd allowed, include more SNPs, or reduce eqtl.n.")
})
