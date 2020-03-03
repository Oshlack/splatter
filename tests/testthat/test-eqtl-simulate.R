context("eQTLSimulate")

data(ex_gff)
data(ex_snps)
eQTLparams <- neweQTLParams()

test_that("eQTLSimulate output is valid", {
    eqtl <- eQTLSimulate()
    expect_true(validObject(eqtl))
    expect_false(any(is.na(eqtl)))
    expect_false(any(sapply(eqtl, is.infinite)))
})

test_that("input data checks", {
    gff_bad <- ex_gff
    gff_bad$V3 <- 'na'
    expect_error(eQTLSimulate(gff = gff_bad),
                 "GFF file did not contain gene features or other issue with 
            file format. See example data.")
    
    vcf_bad <- ex_snps
    vcf_bad$V1 <- NULL
    expect_error(eQTLsnps(vcf_bad, eQTLparams),
                 "snps not in the expected VCF format. See example data.")
})


