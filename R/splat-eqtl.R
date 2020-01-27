#' Splat-eQTL
#'
#' Simulate mean gene counts for a population of samples, such that a defined 
#' number of associations (i.e. cis-eQTL) between markers (i.e. eSNPs) and genes 
#' (i.e. eGenes) exist.
#' 
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details.
#' @param gff_file Path to a GFF/GFT file containing genes to include.
#' @param snp_file Path to real/simulated genotype data in .vcf format. 
#'         Where each column is a sample and each row is a SNP.
#' @param eqtlES_shape Effect Size shape parameter (default estimated from 
#' GTEx thyroid cis-eQTL data)
#' @param eqtlES_rate Effect Size rate parameter (default estimated from 
#' GTEx thyroid cis-eQTL data)
#' @param esnp.n Number of eSNP-eQTL associations to include
#' @param eqtl.dist Distance between eSNP and eGene (TSS). 
#' @param eqtl.maf Desired Minor Allele Frequency (MAF) of eSNPs to include
#' @param eqtl.mafd Maximum variation from eqtl.maf to include as eSNP
#' @param eqtl.bcv Biological Coefficient of Variation (default = 0.4)
#' @param eqtl.save logical. Whether to save eQTL key and mean matrix.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'        
#' @details
#' Parameters can be set in a variety of ways. If no parameters are provided
#' the default parameters are used. Any parameters in \code{params} can be
#' overridden by supplying additional arguments through a call to
#' \code{\link{setParams}}. This design allows the user flexibility in
#' how they supply parameters and allows small adjustments without creating a
#' new \code{SplatParams} object. See examples for a demonstration of how this
#' can be used.
#' 
#' The eQTL Gene Mean simulation involves the following steps:
#' \enumerate{
#'     \item Load and format gene (GFF/GTF) and SNP (genotype) data.
#'     \item Select eGenes-eSNPs pairs and assign effect sizes.
#'     \item Generate normalized gene mean expression matrix for the population. 
#'     \item Set a gene mean expression value (not normalized) for each gene. 
#'     \item Generate a gene mean expression matrix for the population.
#'     \item (optional) Save eQTL key (i.e. pairs)
#' }
#'
#' @return GeneMeansPop Matrix containing the simulated mean gene expression
#' value for each gene (row) and each sample in the population (column).
#' intermediate values.
#' 
#' @seealso
#' \code{\link{splateQTLgenes}}, \code{\link{splateQTLsnps}},
#' \code{\link{splateQTLpairs}}, \code{\link{splateQTLnormMeansMatrix}},
#' \code{\link{splateQTLGeneMeans}}, \code{\link{splateQTLMeansMatrix}}

#' @examples
#' # Load example data
#' library(scater)
#' set.seed(1)
#' sce <- mockSCE()
#' params <- splatEstimate(sce)
#' pop.gMeans <- splateQTL(params)
#' 
#' @export
#' @importFrom utils write.table
splateQTL <- function(params = newSplatParams(),
                      gff_file = 'test_data/test.gff3', 
                      snp_file = 'test_data/test.vcf', 
                      eqtlES_shape = 2.740558,  # Estimated from GTEx thyroid cis-eQTL data
                      eqtlES_rate = 6.441281,  # Estimated from GTEx thyroid cis-eQTL data 
                      esnp.n = 100, 
                      eqtl.dist = 1000000, 
                      eqtl.maf = 0.1, 
                      eqtl.mafd = 0.01,
                      eqtl.bcv = 0.4, # Estimated from GTEx thyroid gene count data
                      eqtl.save = TRUE, ...) {
    
    # Load and format gene (GFF/GTF) and SNP (genotype) data.
    genes <- splateQTLgenes(gff_file)
    snps <- splateQTLsnps(snp_file, eqtl.maf, eqtl.mafd)

    # Select eGenes-eSNPs pairs and assign effect sizes.
    pairs <- splateQTLpairs(genes, snps, esnp.n, eqtl.dist, eqtlES_shape, eqtlES_rate)
    
    # Generate normalized gene mean expression matrix for the population.
    nGeneMeansPop <- splateQTLnormMeansMatrix(snps, pairs)
    
    # Set a gene mean expression value (not normalized) for each gene. 
    pairs <- splateQTLGeneMeans(params, pairs, eqtl.bcv)
    
    # Generate a gene mean expression matrix for the population.
    GeneMeansPop <- splateQTLMeansMatrix(pairs, nGeneMeansPop)
    
    # (optional) Save eQTL key (i.e. pairs)
    if (eqtl.save) {
        dir.create('eqtl_out', showWarnings = FALSE)
        now <- Sys.time()
        write.table(pairs, paste0("eqtl_out/", format(now, "%Y%m%d_%H%M_"), 
                                  "eQTL_key.csv"),
                    sep=' ', quote = FALSE, row.names = FALSE)
    }
    
    return(GeneMeansPop)
}

#' Process gene data
#' 
#' Read in GFF/GTF file (ignorning header) and select only sequences where the 
#' feature is listed as a gene. Then get the Transcriptional Start Site for 
#' each gene (depending on strand direction).
#'
#' @param gff_file Path to GFF/GTF file
#' 
#' @return A dataframe containing gene IDs and locations.
#' @importFrom utils read.delim
splateQTLgenes <- function(gff_file){
    
    gff <- read.delim(gff_file, header=F, comment.char='#')
    genes <- gff[gff[,3]=="gene",]
    genes$TSS <- genes[,4]  #  + strand genes
    genes$TSS[genes[,7] == '-'] <- genes[,5][genes[,7] == '-'] #  - strand genes
    genes$geneID <- c(paste0("gene", 1:dim(genes)[1]))
    genes <- genes[,c('geneID','TSS')]
    
    return(genes)
}


#' Process genotype data
#' 
#' Read in SNP (genotype matrix) file and remove extra columns. Then calculate
#' the Minor Allele Frequency and filter SNPs outside of the MAF range defined.
#'
#' @param snp_file Path to real/simulated genotype data in .ped.gen format
#' @param eqtl.maf Desired Minor Allele Frequency (MAF) of eSNPs to include
#' @param eqtl.mafd Maximum variation from eqtl.maf to include as eSNP
#'
#' @return A dataframe containing SNP names, locations, and sample genotypes.
#' @importFrom utils read.delim
splateQTLsnps <- function(snp_file, eqtl.maf, eqtl.mafd){
    
    # Read in genotype matrix from HAPGEN2 and reformat 
    snps <- read.delim(snp_file, header=F, comment.char='#')
    snps[, c('V1', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9')] <- NULL
    names(snps) <- c('loc', paste0('Sample', 1:(dim(snps)[2]-1)))
    snps[] <- lapply(snps, function(x) gsub("0/0", 0.0, x))
    snps[] <- lapply(snps, function(x) gsub("0/1", 1.0, x))
    snps[] <- lapply(snps, function(x) gsub("1/1", 2.0, x))
    snps <- as.data.frame(sapply(snps, as.numeric))
    snps <- cbind(eSNP = paste0('snp', snps$loc), snps, stringsAsFactors=FALSE)
    
    # Filter out SNPs not within MAF requested
    snps$MAF <- rowSums(snps[,3:dim(snps)[2]]) / ((dim(snps)[2]-2) * 2)
    snps <- subset(snps, MAF > eqtl.maf-eqtl.mafd &
                             MAF < eqtl.maf+eqtl.mafd)
    
    return(snps)
}


#' Select eGenes-eSNPs pairs and assign effect sizes.
#' 
#' Randomly pairs N eSNPs with an eGene within the designated window size 
#' (eqtl.dist) and assigns each pair an effect size sampled from a gamma 
#' distribution parameterized using the effect sizes from a builk eQTL study 
#' using the GTEx data from the thyroid tissue.
#'
#' @param genes Dataframe with gene ID and location
#' @param snps Dataframe with SNP ID, location, and sample genotypes
#' @param eqtlES_shape Effect Size shape parameter (default estimated from 
#' GTEx thyroid cis-eQTL data)
#' @param eqtlES_rate Effect Size rate parameter (default estimated from 
#' GTEx thyroid cis-eQTL data)
#' @param esnp.n Number of eSNP-eQTL associations to include
#' @param eqtl.dist Distance between eSNP and eGene (TSS). 
#' 
#' @return A dataframe eSNPs-eGenes pair assignments and their effect sizes 
#'
splateQTLpairs <- function(genes, snps, esnp.n, eqtl.dist, eqtlES_shape, 
                           eqtlES_rate){
    
    # Set up dataframe to save info about selected eSNP-eGENE pairs 
    snps_list <- snps$eSNP
    pairs <- data.frame(genes)
    pairs$eSNP <- NaN
    pairs$EffectSize <- 0
    genes2 <- data.frame(genes)
    
    for(i in 1:esnp.n){
        again <- TRUE
        while (again == TRUE){
            s <- sample(snps_list, 1)
            snps_list <- snps_list[!snps_list==s]
            l <- snps[snps$eSNP == s, ]$loc
            matches <- subset(genes2, TSS > l - eqtl.dist & TSS < l + eqtl.dist)
            if(dim(matches)[1] > 0){
                match <- sample(matches$geneID, 1)
                again <- FALSE
            }
        }
        
        genes2 <- genes2[!(genes2$geneID==match),]
        ES <- rgamma(1, shape = eqtlES_shape, rate = eqtlES_rate)
        
        pairs[pairs$geneID == match, ]$eSNP <- s
        pairs[pairs$geneID == match, ]$EffectSize <- ES
        
        if(length(snps_list) == 0) {
            warning("Could not find n eSNPs within eqtl.dist of genes provided
                    in the gff file. Either increase the eqtl.dist window or 
                    include more SNPs.")
        }
    }
    
    return(pairs)
}


#' Generate normalized mean gene expression matrix for whole eQTL population
#'
#' Use the approach outlined in Huang et. al 2018 (NAR) to assign normalized
#' mean expression levels for each gene for each sample. Where:
#' y = Effect Size * genotype + error, where error ~ Norm(0,1)
#'
#' @param pairs A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' @param snps The dataframe with the genetic marker info
#'
#' @return normGeneMeansPop: matrix of normalized mean gene exp. levels.
#'
#' @importFrom stats rnorm
splateQTLnormMeansMatrix <- function(snps, pairs) {
    
    # Generate matrix of normalized mean expression values
    samples <- names(snps)[grepl('Sample', names(snps))]
    norm_matrix <- data.frame(matrix(ncol=length(samples), nrow=dim(pairs)[1], 
                                     dimnames=list(pairs$geneID, samples)))
    
    for(g in pairs$geneID) { 
        ES <- pairs[pairs$geneID==g, 'EffectSize']
        
        for(s in samples) {
            error <- rnorm(1, mean = 0, sd = 1)
            if(ES != 0){
                eSNPsample <- pairs[pairs$geneID==g, 'eSNP']
                genotype <- snps[snps$eSNP==eSNPsample,][[s]]
            } else{ 
                genotype <- 0 # just assign 0 because ES = 0 so it won't matter
            }
            norm_matrix[g,s] <- (ES * genotype) + error
        }
    }
    
    return(norm_matrix)
}


#' Set a gene mean expression value (not normalized) for each gene.
#'
#' Assign a mean expression value to each gene, sampled from a gamma 
#' distribution parameterized by splatEstimate, then calculate the BCV and
#' the standard deviation
#'
#' @param params Output from splateQTL_DefineAssociations
#' @param pairs A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' 
#' @return the eSNP-eGene pairs dataframe updated to include the mean gene 
#' expression level, BCV, and standard deviation.
#'
splateQTLGeneMeans <- function(params, pairs, eqtl.bcv){
    
    # Load parameters generated from real data using splatEstimate()
    nGenes <- dim(pairs)[1]
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    out.prob <- getParam(params, "out.prob")
    out.facLoc <- getParam(params, "out.facLoc")
    out.facScale <- getParam(params, "out.facScale")
    #bcv.common <- getParam(params, "bcv.common")
    #bcv.df <- getParam(params, "bcv.df")

    # Simulate base gene means then add outliers
    base.means.gene <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)
    outlier.facs <- getLNormFactors(nGenes, out.prob, 0, out.facLoc,
                                    out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]
    
    # Calculate BCV and thus standard deviation
    pairs$expMean <- means.gene
    #if (is.finite(bcv.df)) {
   #     pairs$BCV <- sqrt((bcv.common + (1 / sqrt(pairs$expMean))) *
     #       (bcv.df / rchisq(1, df = bcv.df)))
    #} else {
   #     warning("'bcv.df' is infinite. This parameter will be ignored.")
    #    pairs$BCV <- (bcv.common + (1 / sqrt(pairs$expMean)))
   # }
    #pairs$Std <- pairs$expMean * pairs$BCV
    pairs$expStd <- pairs$expMean * eqtl.bcv
    
    return(pairs)
}


#' Un-normalize the mean gene expression matrix 
#'
#' @param params A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' @param nMeans The normalized gene expression means for the population
#' 
#' @details
#' For each gene/sample, the normalized expression value (from rnorm) is 
#' transformed to the cumulative density function (pnorm) between 0 and 1, this
#' value is then inversed (qnorm) to map the probability to a value defined by
#' the gene mean assigned in splateQTLGeneMeans.
#'
#' @return MeansMatrix: Matrix of simulated gene means for eQTL population.
#' 
#' @importFrom stats pnorm qnorm
#' 
splateQTLMeansMatrix <- function(pairs, nMeans){
    
    MeansMatrix <- data.frame(nMeans)
    
    for(g in row.names(MeansMatrix)){
        mean.qn <- pairs[pairs$geneID == g,]$expMean
        std.qn <- pairs[pairs$geneID == g,]$expStd
        for(s in names(MeansMatrix)){
            rnorm.tmp <- MeansMatrix[g, s]
            pnorm.tmp <- pnorm(rnorm.tmp, 0, 1)
            MeansMatrix[g, s] <- qnorm(pnorm.tmp, mean.qn, std.qn)
        }
    }
    
    MeansMatrix[MeansMatrix < 0] <- 1e-08
    
    return(MeansMatrix)
}
