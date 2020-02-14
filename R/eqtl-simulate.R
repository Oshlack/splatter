#' eQTL-Simulate
#'
#' Simulate mean gene counts for a population of samples, such that a defined 
#' number of associations (cis-eQTL) between markers (eSNPs) and genes 
#' (eGenes) exist.
#' 
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details.
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' @param gff Dataframe containing the genes to include in GFF/GTF format.
#' @param vcf Dataframe containing real/simulated genotype data in .vcf format. 
#'         Where each column is a sample and each row is a SNP.
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
#'     \item (optional) Save eQTL key (pairs)
#' }
#'
#' @return GeneMeansPop Matrix containing the simulated mean gene expression
#' value for each gene (row) and each sample in the population (column).
#' intermediate values.
#' 
#' @seealso
#' \code{\link{eQTLgenes}}, \code{\link{eQTLsnps}},
#' \code{\link{eQTLpairs}}, \code{\link{eQTLnormMeansMatrix}},
#' \code{\link{eQTLGeneMeans}}, \code{\link{eQTLMeansMatrix}}

#' @examples
#' # Load example data
#' library(scater)
#' set.seed(1)
#' sce <- mockSCE()
#' params <- splatEstimate(sce)
#' eqtlparams <- eQTLEstimate()
#' data(ex_gff)
#' data(ex_snps)
#' pop.gMeans <- eQTLSimulate(params=params, eQTLparams=eqtlparams)
#' 
#' @export
#' @importFrom utils write.table
eQTLSimulate <- function(params = newSplatParams(),
                      eQTLparams = neweQTLParams(),
                      gff = ex_gff, 
                      vcf = ex_snps, 
                      eqtl.save = TRUE, ...) {
    
    # Load and format gene (GFF/GTF) and SNP (genotype) data.
    genes <- eQTLgenes(gff)
    snps <- eQTLsnps(vcf, eQTLparams)

    # Select eGenes-eSNPs pairs and assign effect sizes.
    pairs <- eQTLpairs(genes, snps, eQTLparams)
    
    # Generate normalized gene mean expression matrix for the population.
    nGeneMeansPop <- eQTLnormMeansMatrix(snps, pairs)
    
    # Select a gene mean expression value & variance for each gene. 
    pairs <- eQTLGeneMeans(params, pairs, eQTLparams)
    
    # Generate a gene mean expression matrix for the population.
    GeneMeansPop <- eQTLMeansMatrix(pairs, nGeneMeansPop, params)
    
    # (optional) Save eQTL key (pairs)
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
#' Read in GFF/GTF file (ignoring header) and select only sequences where the 
#' feature is listed as a gene. Then get the Transcriptional Start Site for 
#' each gene (depending on strand direction).
#'
#' @param gff Dataframe containing the genes to include in GFF/GTF format.
#' 
#' @return A dataframe containing gene IDs and locations.
#' @importFrom utils read.delim
eQTLgenes <- function(gff){
    
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
#' @param vcf Dataframe containing real/simulated genotype data in .vcf format. 
#'         Where each column is a sample and each row is a SNP.
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#'
#' @return A dataframe containing SNP names, locations, and sample genotypes.
#' @importFrom utils read.delim
eQTLsnps <- function(vcf, eQTLparams){
    
    eqtl.maf <- getParam(eQTLparams, "eqtl.maf")
    eqtl.mafd <- getParam(eQTLparams, "eqtl.mafd")
    MAF <- NULL  # locally binding variables
    
    # Read in genotype matrix in .vcf format
    #snps <- read.delim(snp_file, header=F, comment.char='#')
    vcf[, c('V1', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9')] <- NULL
    names(vcf) <- c('loc', paste0('Sample', 1:(dim(vcf)[2]-1)))
    vcf[] <- lapply(vcf, function(x) gsub("0/0", 0.0, x))
    vcf[] <- lapply(vcf, function(x) gsub("0/1", 1.0, x))
    vcf[] <- lapply(vcf, function(x) gsub("1/1", 2.0, x))
    vcf <- as.data.frame(sapply(vcf, as.numeric))
    vcf <- cbind(eSNP = paste0('snp', vcf$loc), vcf, stringsAsFactors=FALSE)
    
    # Filter out SNPs not within MAF requested
    vcf$MAF <- rowSums(vcf[,3:dim(vcf)[2]]) / ((dim(vcf)[2]-2) * 2)
    snps <- subset(vcf, MAF > eqtl.maf-eqtl.mafd &
                             MAF < eqtl.maf+eqtl.mafd)
    if(dim(snps)[1] < getParam(eQTLparams, "eqtl.n")){
        warning("Not enough SNPs within desired MAF range. Either increase the
                    eqtl.mafd allowed, include more SNPs, or reduce eqtl.n.")
    }
    
    return(snps)
}


#' Select eGenes-eSNPs pairs and assign effect sizes.
#' 
#' Randomly pairs N eSNPs with an eGene within the designated window size 
#' (eqtl.dist) and assigns each pair an effect size sampled from a gamma 
#' distribution parameterized using the effect sizes from a bulk eQTL study 
#' using the GTEx data from the thyroid tissue.
#'
#' @param genes Dataframe with gene ID and location
#' @param snps Dataframe with SNP ID, location, and sample genotypes
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' 
#' @return A dataframe eSNPs-eGenes pair assignments and their effect sizes 
#'
eQTLpairs <- function(genes, snps, eQTLparams){
    eqtl.n <- getParam(eQTLparams, "eqtl.n")
    eqtl.dist <- getParam(eQTLparams, "eqtl.dist")
    eqtlES_shape <- getParam(eQTLparams, "eqtlES.shape") 
    eqtlES_rate <- getParam(eQTLparams, "eqtlES.rate")
    TSS <- NULL  # locally binding variables
    
    # Set up dataframe to save info about selected eSNP-eGENE pairs 
    snps_list <- snps$eSNP
    pairs <- data.frame(genes)
    pairs$eSNP <- NaN
    pairs$EffectSize <- 0
    genes2 <- data.frame(genes)
    
    for(i in 1:eqtl.n){
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
                    in the GFF file. Either increase the eqtl.dist window or 
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
#' @param snps The dataframe with the genetic marker info
#' @param pairs A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' 
#' @return normGeneMeansPop: matrix of normalized mean gene exp. levels.
#'
#' @importFrom stats rnorm
eQTLnormMeansMatrix <- function(snps, pairs) {
    
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
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details.
#' @param pairs A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' 
#' @return the eSNP-eGene pairs dataframe updated to include the mean gene 
#' expression level, BCV, and standard deviation.
#' @importFrom data.table between
eQTLGeneMeans <- function(params, pairs, eQTLparams){
    
    # Load parameters generated from real data using splatEstimate()
    nGenes <- dim(pairs)[1]
    out.prob <- getParam(params, "out.prob")
    out.facLoc <- getParam(params, "out.facLoc")
    out.facScale <- getParam(params, "out.facScale")
    bulk_mean_shape <- getParam(eQTLparams, "bulkmean.shape")
    bulk_mean_rate <- getParam(eQTLparams, "bulkmean.rate")
    cv.param <- getParam(eQTLparams, "bulkcv.param")

    # Sample gene means and add outliers
    base.means.gene <- rgamma(nGenes, shape = bulk_mean_shape, 
                              rate = bulk_mean_rate)
    outlier.facs <- getLNormFactors(nGenes, out.prob, 0, out.facLoc,
                                    out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]
    pairs$expMean <- means.gene
    pairs$expCV <- NULL
    
    # Sample coefficient of variation for each gene
    for (g in 1:nrow(pairs)){
        mean <- pairs[g, 'expMean']
        bin <- cv.param[(cv.param$start < mean) & (cv.param$end >= mean), ]
        cv_shape <- bin$shape
        cv_rate <- bin$rate
        pairs[g,'expCV'] <- rgamma(1, shape = cv_shape, rate = cv_rate)
    }

    return(pairs)
}


#' Project normalized gene expression matrix into mean gene expression matrix 
#'
#' @param pairs A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' @param nGeneMeansPop The normalized gene expression means for the population
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details.
#' 
#' @details
#' For each gene/sample, the normalized expression value (from rnorm) is 
#' transformed to the cumulative density function (pnorm) between 0 and 1, this
#' value is then quantile normalized (qgamma) using the gamma distribution 
#' parameterized from splatEstimate(). 
#'
#' @return MeansMatrix: Matrix of simulated gene means for eQTL population.
#' 
#' @importFrom stats pnorm qnorm qgamma
#' 
eQTLMeansMatrix <- function(pairs, nGeneMeansPop, params){
    
    # For each gene, simulate normal dist. of mean expression across population.
    MeansMatrix <- data.frame(nGeneMeansPop)
    for(g in row.names(MeansMatrix)){
        mean.gene <- pairs[pairs$geneID == g,]$expMean
        cv.gene <- pairs[pairs$geneID == g,]$expCV
        sd.gene <- cv.gene * mean.gene
        norm.mean <- mean(unlist(MeansMatrix[g, ]))
        norm.sd <- sd(unlist(MeansMatrix[g, ]))
        for(s in names(MeansMatrix)){
            rnorm.tmp <- MeansMatrix[g, s]
            pnorm.tmp <- pnorm(rnorm.tmp, norm.mean, norm.sd)
            MeansMatrix[g, s] <- qnorm(pnorm.tmp, mean.gene, sd.gene)
        }
    }
    MeansMatrix[MeansMatrix < 0] <- 0
   
     # For each sample, quantile gamma normalize expression across genes.
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    for(s in names(MeansMatrix)){
        sample_mean <- mean(MeansMatrix[, s])
        sample_sd <- sd(MeansMatrix[, s])
        pnorm.tmp <- pnorm(MeansMatrix[, s], sample_mean, sample_sd)
        MeansMatrix[, s] <- qgamma(pnorm.tmp, shape=mean.shape,
                                  rate=mean.rate)
    }
    
    return(MeansMatrix)
}


#utils::globalVariables(c("ex_gff", "ex_snps"))