context("eQTLSimulate")

data(ex_gff)
data(ex_snps)
vcf <- ex_snps[1:5000,]
eparams <- setParams(neweQTLParams(), list(eqtl.n = 10, seed = 42))

test_that("eQTLSimulate output is valid and works", {
    eqtl <- eQTLSimulate(eQTLparams = eparams, 
                         vcf = vcf, eqtl.save = FALSE)
    expect_true(validObject(eqtl))
    expect_false(any(is.na(eqtl)))
    expect_false(any(sapply(eqtl, is.infinite)))
    expect_length(eqtl, 100)
    #expect_true(round(eqtl[1,1], 2) == 30.68)
})

test_that("input data checks", {
    gff_bad <- ex_gff
    gff_bad$V3 <- 'na'
    expect_error(eQTLSimulate(gff = gff_bad),
                 "GFF file did not contain gene features or other issue with 
            file format. See example data.")
    
    vcf_bad <- ex_snps
    vcf_bad$V1 <- NULL
    expect_error(eQTLsnps(vcf_bad, eparams),
                 "snps not in the expected VCF format. See example data.")
    
    vcf_bad2 <- ex_snps[1:5,]
    expect_error(eQTLSimulate(vcf = vcf_bad2), 
    "Not enough SNPs within desired MAF range. Increase the
                    eqtl.mafd allowed, include more SNPs, or reduce eqtl.n.")
})
