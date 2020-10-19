#' Generate mock gff 
#'
#' Quick function to generate a mock gff. 
#' 
#' @param n.genes Number of genes in mock gff file
#' @param chromosome Chromosome name 
#' 
#' @return data.frame containing mock gff data.
#' 
#' @export
#' 
mockGFF <- function(n.genes = 500, chromosome = 22){
    
    mock.gff <- data.frame(list(V1 = chromosome,
                                V2 = "source",
                                V3 = "gene",
                                V4 = sort(sample(1e4:2e8, n.genes))))
    mock.gff$V5 <- mock.gff$V4 + floor(rnorm(n.genes, 1500, 1000))
    mock.gff[, c("V6", "V7", "V8", "V9")] <- "."
    
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
#' 
#' @return data.frame containing mock gff data.
#'
#' @importFrom stats setNames
#' 
#' @export
#' 
mockVCF <- function(n.snps = 1e4, n.samples = 10, chromosome = 22){
    
    
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
        stop("Creating a mock VCF requires the 'VariantAnnotation' package.")
    }
    
    sample_names <- paste0("sample_", formatC(1:n.samples, 
                                              width = nchar(n.samples),
                                              format = "d", 
                                              flag = "0"))
    snp_names <- paste0("snp_", formatC(1:n.snps, 
                                        width = nchar(n.snps),
                                        format = "d",
                                        flag = "0"))
    # rowRanges
    vcf.rowRanges <- GenomicRanges::GRanges(
        seqnames = Rle(rep(chromosome, n.snps)),
        ranges = IRanges(sample(1:2e8, n.snps, replace = FALSE), 
                         names = snp_names),
        strand = Rle(strand(rep("*", n.snps))),
        paramRangeID = Rle(rep(NA, n.snps)))
    
    # Geno
    genotypes <- c(rep("0|0", 100), rep("1|0", 10), rep("1|1", 10))
    geno <- setNames(data.frame(matrix(ncol = n.samples, nrow = n.snps)), 
                     sample_names)
    geno[,  sample_names] <- sample(genotypes, n.samples * n.snps, 
                                    replace = TRUE)
    row.names(geno) <- snp_names
    geno <- as.matrix(geno)
    
    # colData, info, and fixed
    vcf.col.data <- DataFrame(list(Samples = 1:n.samples))
    row.names(vcf.col.data) <- sample_names
    
    vcf.info <- DataFrame(list(VT = rep("SNP", n.snps)))
    row.names(vcf.info) <- snp_names
    
    vcf.fixed <- DataFrame(
        REF = Biostrings::DNAStringSet(rep("A", n.snps)),
        ALT = rep(Biostrings::DNAStringSetList("T"), n.snps),
        QUAL = rep(100, n.snps),
        FILTER = rep("PASS", n.snps))
    
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
#' 
#' @return matrix containing mock bulk expression data.
#' 
#' @export

mockBulkMatrix <- function(n.genes = 1000, n.samples = 100){
    
    tmp.params <- newSplatPopParams()
    mean.shape <- getParam(tmp.params, "pop.mean.shape")
    mean.rate <- getParam(tmp.params, "pop.mean.rate")
    cv.df <- getParam(tmp.params, "pop.cv.param")
    cv.shape <- cv.df[5, "shape"]
    cv.rate <- cv.df[5, "rate"]
    
    key <- data.frame(list(id = c(1:n.genes),
                           mean = rgamma(n.genes, mean.shape, mean.rate),
                           cv = rgamma(n.genes, cv.shape, cv.rate)))
    
    mock.means <- lapply(key$id, function(g) rnorm(n.samples,
                                              mean = key[key$id == g,]$mean,
                                              sd = key[key$id == g,]$mean *
                                                  key[key$id == g,]$cv))
    
    mock.means <- do.call(rbind, mock.means)
    mock.means[mock.means < 0] <- 0
    
    return(mock.means)
}

#' Generate mock eQTL mapping results
#'
#' Quick function to generate mock eQTL mapping results, with parameters 
#' estimated using real eQTL mapping results from GTEX using thyroid tissue. 
#' 
#' @param n.genes Number of genes in mock eQTL data.
#' 
#' @return data.frame containing mock bulk eQTL mapping results.
#' 
#' @export
#' 
mockBulkeQTL <- function(n.genes = 1000){
    
    tmp.params <- newSplatPopParams()
    eqtl.shape <- getParam(tmp.params, "eqtl.ES.shape")
    eqtl.rate <- getParam(tmp.params, "eqtl.ES.rate")
    
    mock.eq <- data.frame(list(gene_id = 1:n.genes, 
                               pval_nominal = 0.01,
                               slope = rgamma(n.genes, eqtl.shape, eqtl.rate)))
    
    return(mock.eq)
}
