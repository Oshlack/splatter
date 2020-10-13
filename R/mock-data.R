#' Generate mock gff 
#'
#' Quick function to generate a mock gff. 
#' 
#' @param n.genes Number of genes in mock gff file
#' @param chromosome Chromosome name 
#' 
#' @return data.frame containing mock gff data from chromosome 1.
#' 
#' @export
#' 
mockGFF <- function(n.genes = 1000, chromosome = 1){
    
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
#' @return data.frame containing mock gff data from chromosome 1.
#'
#' @export
#' 
mockVCF <- function(n.snps = 1e4, n.samples = 10, chromosome = 1){
    
    if (!requireNamespace("vcfR", quietly = TRUE)) {
        stop("Creating a mock VCF requires the 'vcfR' package.")
    }
    data(vcfR_example, envir = environment())
    vcf@meta[2] <- "##source=\"Mock\""    
    
    fix <- t(do.call(rbind, list(CHROM = rep(chromosome, n.snps),
                              POS = sort(sample(1:2e8, n.snps, 
                                                replace = FALSE)),
                              ID = NA,
                              REF = "A",
                              ALT = "T", 
                              QUAL = "50.0",
                              FILTER = NA,
                              INFO = NA)))
    
    genotypes <- as.data.frame(vcf@gt)
    format <- genotypes$FORMAT[1]
    genotypes$FORMAT <- NULL
    genotypes <- na.omit(as.vector(t(genotypes)))
    
    gt <- data.frame(list(FORMAT = rep(format, n.snps)))
    
    samples <- paste0("sample_", formatC(1:n.samples, 
                                         width = nchar(n.samples),
                                         format = "d", flag = "0"))
    
    gt[,  samples] <- sample(genotypes, n.samples * n.snps, replace = TRUE)
    gt <- as.matrix(gt)
    
    vcf@fix <- fix
    vcf@gt <- gt

    return(vcf)
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
