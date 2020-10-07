#' Generate mock gff 
#'
#' Quick function to generate a mock gff. 
#' 
#' @param n_genes Number of genes in mock gff file
#' @param chromosome Chromosome value. Default=1.
#' 
#' @return data.frame containing mock gff data from chromosome 1.
#' 
#' @export
#' 
mock_gff <- function(n_genes=1000, chromosome=1){
    
    mock_gff <- data.frame(list(V1 = chromosome,
                                V2 = "source",
                                V3 = "gene",
                                V4 = sort(sample(1e4:2e8, n_genes))))
    mock_gff$V5 <- mock_gff$V4 + floor(rnorm(n_genes, 1500, 1000))
    mock_gff[, c("V6", "V7", "V8", "V9")] <- "."
    
    return(mock_gff)
}


#' Generate mock vcf 
#'
#' Quick function to generate mock vcf file. Note this data has unrealistic
#' population structure.
#' 
#' @param n_snps Number of SNPs in mock vcf file. Default=10,000.
#' @param n_samples Number of samples in mock bulk data. Default=10.
#' @param chromosome Chromosome value. Default=1.
#' 
#' @return data.frame containing mock gff data from chromosome 1.
#'
#' @export
#' 
mock_vcf <- function(n_snps=1e4, n_samples=10, chromosome=1){
    
    if (requireNamespace("vcfR", quietly = TRUE))
    data(vcfR_example)
    mock_vcf <- vcf
    rm(gff, vcf, dna, pos=1)
    mock_vcf@meta[2] <- "##source=\"Mock\""    
    
    fix <- t(do.call(rbind, list(CHROM = rep(chromosome, n_snps),
                              POS = sort(sample(1:2e8, n_snps, replace=FALSE)),
                              ID = NA,
                              REF = "A",
                              ALT = "T", 
                              QUAL = "50.0",
                              FILTER = NA,
                              INFO = NA)))
    
    genotypes <- as.data.frame(mock_vcf@gt)
    format <- genotypes$FORMAT[1]
    genotypes$FORMAT <- NULL
    genotypes <- na.omit(as.vector(t(genotypes)))
    
    gt <- data.frame(list(FORMAT = rep(format, n_snps)))
    
    samples <- paste0("sample_", formatC(1:n_samples, 
                                         width= nchar(n_samples),
                                         format="d", flag="0"))
    
    gt[,  samples] <- sample(genotypes, n_samples * n_snps, replace=TRUE)
    gt <- as.matrix(gt)
    
    mock_vcf@fix <- fix
    mock_vcf@gt <- gt

    return(mock_vcf)
}


#' Generate mock bulk population scale expression data
#'
#' Quick function to generate mock bulk expression data for a population, with
#' parameters estimated using real thyroid tissue data from GTEx. 
#' 
#' @param n_genes Number of genes in mock bulk data. Default=1,000.
#' @param n_samples Number of samples in mock bulk data. Default=100.
#' 
#' @return data.frame containing mock bulk expression data.
#' 
#' @export

mock_bulk_matrix <- function(n_genes=1000, n_samples=100){
    
    tmp_params <- newSplatPopParams()
    mean_shape <- getParam(tmp_params, "pop.mean.shape")
    mean_rate <- getParam(tmp_params, "pop.mean.rate")
    cv_df <- getParam(tmp_params, "pop.cv.param")
    cv_shape <- cv_df[5, "shape"]
    cv_rate <- cv_df[5, "rate"]
    
    key <- data.frame(list(id = c(1:n_genes),
                           mean = rgamma(n_genes, mean_shape, mean_rate),
                           cv = rgamma(n_genes, cv_shape, cv_rate)))
    
    mock_means <- lapply(key$id, function(g) rnorm(n_samples,
                                              mean = key[key$id == g,]$mean,
                                              sd = key[key$id == g,]$mean *
                                                  key[key$id == g,]$cv))
    
    mock_means <- data.frame(do.call(rbind, mock_means))
    mock_means[mock_means < 0] <- 0
    
    return(mock_means)
}

#' Generate mock eQTL mapping results
#'
#' Quick function to generate mock eQTL mapping results, with parameters 
#' estimated using real eQTL mapping results from GTEX using thyroid tissue. 
#' 
#' @param n_genes Number of genes in mock eQTL data. Default=1,000.
#' 
#' @return data.frame containing mock bulk eQTL mapping results.
#' 
#' @export
#' 
mock_bulk_eqtl <- function(n_genes=1000){
    
    tmp_params <- newSplatPopParams()
    eqtl_shape <- getParam(tmp_params, "eqtl.ES.shape")
    eqtl_rate <- getParam(tmp_params, "eqtl.ES.rate")
    
    mock_eq <- data.frame(list(gene_id = 1:n_genes, 
                               pval_nominal = 0.01,
                               slope = rgamma(n_genes, eqtl_shape, eqtl_rate)))
    
    return(mock_eq)
}
