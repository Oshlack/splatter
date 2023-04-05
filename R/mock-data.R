#' Generate mock gff
#'
#' Quick function to generate a mock gff.
#'
#' @param n.genes Number of genes in mock gff file
#' @param chromosome Chromosome name
#' @param chr.length Length of mock chromosome
#' @param seed Optional: seed for random seed
#'
#' @return data.frame containing mock gff data.
#'
#' @examples
#' gff <- mockGFF()
#'
#' @export
mockGFF <- function(n.genes = 50, chromosome = 1, chr.length = 2e6, seed = NULL){

    if (is.null(seed)) {
        V4 = sort(sample(1e4:chr.length, n.genes))
        V5 <- V4 + floor(rnorm(n.genes, 1500, 1000))
    } else {
        withr::with_seed(seed, {
            V4 = sort(sample(1e4:chr.length, n.genes))
            V5 <- V4 + floor(rnorm(n.genes, 1500, 1000))
        })
    }

    mock.gff <- data.frame(
        V1 = chromosome,
        V2 = "source",
        V3 = "gene",
        V4 = V4,
        V5 = V5,
        V6 = ".",
        V7 = ".",
        V8 = ".",
        V9 = "."
    )

    return(mock.gff)
}


#' Generate mock vcf
#'
#' Quick function to generate mock vcf file. Note this data has unrealistic
#' population structure.
#'
#' @param n.snps Number of SNPs in mock vcf file.
#' @param n.samples Number of samples in mock bulk data.
#' @param chromosome Chromosome name
#' @param chr.length Length of mock chromosome
#' @param seed Optional: seed for random seed
#'
#' @return data.frame containing mock vcf data.
#'
#' @examples
#' vcf <- mockVCF()
#'
#' @export
#' @importFrom stats setNames
mockVCF <- function(n.snps = 200, n.samples = 5,
                    chromosome = 1,  chr.length = 2e6, seed = NULL){

    checkDependencies(deps = "VariantAnnotation")

    sample_names <- paste0("sample_", formatC(seq_len(n.samples),
                                              width = nchar(n.samples),
                                              format = "d",
                                              flag = "0"))
    snp_names <- paste0("snp_", formatC(seq_len(n.snps),
                                        width = nchar(n.snps),
                                        format = "d", flag = "0"))

    if (is.null(seed)) {
        ranges <- IRanges::IRanges(sample(seq_len(chr.length), n.snps,
                                          replace = FALSE),
                                   names = snp_names)
    } else {
        withr::with_seed(seed, {
            ranges <- IRanges::IRanges(sample(seq_len(chr.length), n.snps,
                                              replace = FALSE),
                                       names = snp_names)
        })
    }

    # rowRanges
    vcf.rowRanges <- GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(rep(chromosome, n.snps)),
        ranges = ranges,
        strand = S4Vectors::Rle(BiocGenerics::strand(rep("*", n.snps))),
        paramRangeID = S4Vectors::Rle(rep(NA, n.snps))
    )

    # Geno
    genotypes <- c(rep("0|0", 100), rep("1|0", 10), rep("1|1", 10))
    geno <- setNames(data.frame(matrix(ncol = n.samples, nrow = n.snps)),
                     sample_names)

    if (is.null(seed)) {
        geno[,  sample_names] <- sample(genotypes, n.samples * n.snps,
                                        replace = TRUE)
    } else {
        withr::with_seed(seed, {
            geno[,  sample_names] <- sample(genotypes, n.samples * n.snps,
                                            replace = TRUE)
        })
    }

    row.names(geno) <- snp_names
    geno <- as.matrix(geno)

    # colData, info, and fixed
    vcf.col.data <- S4Vectors::DataFrame(list(Samples = seq_len(n.samples)))
    row.names(vcf.col.data) <- sample_names

    vcf.info <- S4Vectors::DataFrame(list(VT = rep("SNP", n.snps)))
    row.names(vcf.info) <- snp_names

    vcf.fixed <- S4Vectors::DataFrame(
        REF = Biostrings::DNAStringSet(rep("A", n.snps)),
        ALT = rep(Biostrings::DNAStringSetList("T"), n.snps),
        QUAL = rep(100, n.snps),
        FILTER = rep("PASS", n.snps)
    )

    # Merge into mock VCF object
    mock.vcf <- VariantAnnotation::VCF(rowRanges = vcf.rowRanges,
                colData = vcf.col.data,
                info = vcf.info,
                fixed = vcf.fixed,
                geno = list(GT=geno),
                collapsed = TRUE)

    return(mock.vcf)
}


#' Generate mock bulk population scale expression data
#'
#' Quick function to generate mock bulk expression data for a population, with
#' parameters estimated using real thyroid tissue data from GTEx.
#'
#' @param n.genes Number of genes in mock bulk data.
#' @param n.samples Number of samples in mock bulk data.
#' @param seed Optional: seed for random seed
#'
#' @return matrix containing mock bulk expression data.
#'
#' @examples
#' bulk <- mockBulkMatrix
#'
#' @export
mockBulkMatrix <- function(n.genes = 100, n.samples = 50, seed = NULL){

    tmp.params <- newSplatPopParams()
    mean.shape <- getParam(tmp.params, "pop.mean.shape")
    mean.rate <- getParam(tmp.params, "pop.mean.rate")
    cv.df <- getParam(tmp.params, "pop.cv.param")
    cv.shape <- cv.df[5, "shape"]
    cv.rate <- cv.df[5, "rate"]

    if (is.null(seed)) {
        key <- data.frame(
            id = seq_len(n.genes),
            mean = rgamma(n.genes, mean.shape, mean.rate),
            cv = rgamma(n.genes, cv.shape, cv.rate)
        )

        mock.means <- lapply(key$id, function(g) {
            rnorm(n.samples, mean = key[key$id == g,]$mean,
                  sd = key[key$id == g,]$mean * key[key$id == g,]$cv)
        })
    } else {
        withr::with_seed(seed, {
            key <- data.frame(
                id = seq_len(n.genes),
                mean = rgamma(n.genes, mean.shape, mean.rate),
                cv = rgamma(n.genes, cv.shape, cv.rate)
            )

            mock.means <- lapply(key$id, function(g) {
                rnorm(n.samples, mean = key[key$id == g,]$mean,
                      sd = key[key$id == g,]$mean * key[key$id == g,]$cv)
            })
        })
    }

    mock.means <- do.call(rbind, mock.means)
    mock.means[mock.means < 0] <- 0

    return(mock.means)
}

#' Generate mock eQTL mapping results
#'
#' Quick function to generate mock eQTL mapping results, with parameters
#' estimated using real eQTL mapping results from GTEx using thyroid tissue.
#'
#' @param n.genes Number of genes in mock eQTL data.
#' @param seed Optional: seed for random seed
#'
#' @return data.frame containing mock bulk eQTL mapping results.
#'
#' @examples
#' eqtl <- mockBulkeQTL()
#'
#' @export
mockBulkeQTL <- function(n.genes = 500, seed = NULL){

    tmp.params <- newSplatPopParams()
    eqtl.shape <- getParam(tmp.params, "eqtl.ES.shape")
    eqtl.rate <- getParam(tmp.params, "eqtl.ES.rate")

    if (is.null(seed)) {
        mock.eq <- data.frame(
            gene_id = seq_len(n.genes),
            pval_nominal = 0.01,
            slope = rgamma(n.genes, eqtl.shape, eqtl.rate)
        )
    } else {
        withr::with_seed(seed, {
            mock.eq <- data.frame(
                gene_id = seq_len(n.genes),
                pval_nominal = 0.01,
                slope = rgamma(n.genes, eqtl.shape, eqtl.rate)
            )
        })
    }

    return(mock.eq)
}


#' Generate set of "empirical" mock data
#'
#' Quick function to generate matching mock VCF, bulk expression, and eQTL data,
#' useful for running splatPopEmpiricalMeans
#'
#' @param n.genes Number of genes in mock eQTL data.
#' @param n.snps Number of SNPs in mock vcf file.
#' @param n.samples Number of samples in mock bulk data.
#' @param chromosome Chromosome name
#' @param chr.length Length of mock chromosome
#' @param seed Optional: seed for random seed
#'
#' @return list(gff=mockGFF, vcf=mockVCF, means=mockMEANS, eqtl=mockEQTL)
#'
#' @examples
#' empirical <- mockEmpiricalSet()
#'
#' @export
#'
mockEmpiricalSet <- function(n.genes = 20, n.snps = 1000, n.samples = 10,
                    chromosome = 1,  chr.length = 2e6, seed = NULL){

    mockGFF <- mockGFF(n.genes = n.genes, chromosome = chromosome,
                   chr.length = chr.length, seed = seed)
    mockVCF <- mockVCF(n.snps = n.snps, n.samples = n.samples,
                       chromosome = chromosome,  chr.length = chr.length,
                       seed = seed)
    mockMEANS <- mockBulkMatrix(n.genes = n.genes, n.samples = n.samples,
                                seed = seed)
    mockEQTL <- mockBulkeQTL(n.genes = n.genes, seed = seed)

    mockEQTL$geneID <- paste0("gene_", formatC(seq_len(n.genes),
                                               width = nchar(n.genes),
                                               format = "d", flag = "0"))

    row.names(mockMEANS) <- mockEQTL$geneID
    colnames(mockMEANS) <- colnames(mockVCF)

    vcfDF <- data.frame(SummarizedExperiment::rowRanges(mockVCF))
    row.names(vcfDF) <- rownames(mockVCF)
    vcfDF$MAF <- VariantAnnotation::snpSummary(mockVCF)$a1Freq

    if (is.null(seed)) {
        eSNPs <- sample(rownames(mockVCF), n.genes)
    } else {
        withr::with_seed(seed, {
            eSNPs <- sample(rownames(mockVCF), n.genes)
        })
    }

    vcfDF <- vcfDF[eSNPs, ]

    mockEQTL$snpID <- eSNPs
    mockEQTL$snpLOC <- vcfDF$start
    mockEQTL$snpCHR <- vcfDF$seqnames
    mockEQTL$snpMAF <- vcfDF$MAF

    mockEQTL[mockEQTL$snpMAF < 0.05, c("snpID", "snpLOC", "snpCHR")] <- NA

    return(list(gff=mockGFF, vcf=mockVCF, means=mockMEANS, eqtl=mockEQTL))
}
