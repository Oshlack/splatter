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
#' @param save.name String to save eQTL key under. Default: 
#'         eqtl_out/YYMMDD_HHMM_eQTL_key.csv
#' @param verbose logical. Whether to print progress messages.
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
#' data(ex_gff)
#' data(ex_snps)
#' pop.gMeans <- eQTLSimulate()
#' 
#' @export
#' @importFrom utils write.table
eQTLSimulate <- function(params = newSplatParams(),
                      eQTLparams = neweQTLParams(),
                      gff = ex_gff, 
                      vcf = ex_snps, 
                      eqtl.save = TRUE,
                      save.name = 'default',
                      verbose = TRUE, ...) {
    
    
    # Set random seed
    seed <- getParam(eQTLparams, "seed")
    set.seed(seed)
    
    if (verbose) {message("Loading GFF and SNP data...")}
    genes <- eQTLgenes(gff)
    snps <- eQTLsnps(vcf, eQTLparams)
    
    if (verbose) {message("Assigning expression, variance, & eQTL effects...")} 
    genes <- eQTLGeneMeans(params, genes, eQTLparams)
    groups <- paste0('g', seq(1, getParam(eQTLparams, "eqtl.groups")))
    
    # Assign eQTL effects for each group
    for(id in groups){
        if(id == 'g1'){
            pairs <- eQTLpairs(id, genes, snps, eQTLparams)
        }else{
            pairs[, paste0("eSNP_", id)] <- pairs$eSNP_g1
            pairs[, paste0("EffectSize_", id)] <- pairs$EffectSize_g1
        }
    }
    
    if(length(groups) > 1){
        pairs <- groupSpecificEffects(eQTLparams, pairs, groups)
    }
    
    nGeneMeansPop <- GeneMeansPop <- FALSE
    final <- list()
    
    for(id in groups){
        
        if (verbose) {message(paste0("Simulating group ", id, "..."))}
        nGeneMeansPop <- eQTLnormMeansMatrix(id, snps, pairs, nGeneMeansPop)
        GeneMeansPop <- eQTLMeansMatrix(id, pairs, nGeneMeansPop, GeneMeansPop)
        GeneMeansPop <- quantileNormalizeSC(params, GeneMeansPop)

        final[[id]] <- GeneMeansPop
    }
    
    # Save eQTL key
    if (eqtl.save) {
        dir.create('eqtl_out', showWarnings = FALSE)
        if(save.name == 'default'){
            now <- Sys.time()
            save <- paste0("eqtl_out/", format(now, "%Y%m%d_%H%M_"), "eQTL_key.csv")
        }else(
            save <- paste0('eqtl_out/', save.name)
        )
        
        write.table(pairs, save, sep=',', quote = FALSE, row.names = FALSE)
        if (verbose) {message(paste0("Saved key to: ", save))}
    }
    
    if (verbose) {message("Done!")} 
    
    if(length(groups) == 1){
        return(as.data.frame(final[id]))
    }else{
        return(final)
    }
    
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
    
    # Test input gff file
    if ((length(names(gff)) < 8 | nrow(gff[gff[,3]=="gene",]) < 1)) {
        stop("GFF file did not contain gene features or other issue with 
            file format. See example data.")
    }
    
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
#' @importFrom dplyr mutate_all
#' 
eQTLsnps <- function(vcf, eQTLparams){
    
    eqtl.maf <- getParam(eQTLparams, "eqtl.maf")
    eqtl.mafd <- getParam(eQTLparams, "eqtl.mafd")
    MAF <- NULL  # locally binding variables
    
    # Test input VCF file
    if (!('V1' %in% names(vcf))) {
        stop("snps not in the expected VCF format. See example data.")
    }
    
    # Read in genotype matrix in .vcf format
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
    return(snps)
}


#' Set a gene mean expression value (not normalized) for each gene.
#'
#' Assign a mean expression value to each gene, sampled from a gamma 
#' distribution parameterized by splatEstimate, then sample a coefficient of 
#' variation for each gene from a gamma distribution parameterized by 
#' splatEstimate for genes in 10 mean expression bins.
#'
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details.
#' @param genes A dataframe of genes included in the sim
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' 
#' @return the gene dataframe updated to include the mean gene expression 
#' level and coefficient of variation.
#' @importFrom data.table between
eQTLGeneMeans <- function(params, genes, eQTLparams){
    
    # Load parameters generated from real data using splatEstimate()
    nGenes <- dim(genes)[1]
    bulk_mean_shape <- getParam(eQTLparams, "bulkmean.shape")
    bulk_mean_rate <- getParam(eQTLparams, "bulkmean.rate")
    cv.param <- getParam(eQTLparams, "bulkcv.param")
    
    # Sample gene means (note, outliers added by splatSimulate()
    base.means.gene <- rgamma(nGenes, shape = bulk_mean_shape, 
                              rate = bulk_mean_rate)
    means.gene <- base.means.gene 
    genes$expMean <- means.gene
    genes$expCV <- NULL
    
    # Sample coefficient of variation for each gene
    for (g in 1:nrow(genes)){
        mean <- genes[g, 'expMean']
        bin <- cv.param[(cv.param$start < mean) & (cv.param$end >= mean), ]
        cv_shape <- bin$shape
        cv_rate <- bin$rate
        genes[g,'expCV'] <- rgamma(1, shape = cv_shape, rate = cv_rate)
    }
    
    return(genes)
}


#' Select eGenes-eSNPs pairs and assign effect sizes.
#' 
#' Randomly pairs N eSNPs with an eGene within the designated window size 
#' (eqtl.dist) and assigns each pair an effect size sampled from a gamma 
#' distribution parameterized using the effect sizes from a bulk eQTL study 
#' using the GTEx data from the thyroid tissue.
#'
#' @param id group id name
#' @param genes Dataframe with gene ID and location
#' @param snps Dataframe with SNP ID, location, and sample genotypes
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' 
#' @return A dataframe eSNPs-eGenes pair assignments and their effect sizes 
#'
eQTLpairs <- function(id, genes, snps, eQTLparams){
    eqtl.n <- getParam(eQTLparams, "eqtl.n")
    if (eqtl.n > dim(genes)[1]){
        eqtl.n <- dim(genes)[1]
    }
    eqtl.dist <- getParam(eQTLparams, "eqtl.dist")
    eqtlES_shape <- getParam(eQTLparams, "eqtlES.shape") 
    eqtlES_rate <- getParam(eQTLparams, "eqtlES.rate")
    TSS <- NULL  # locally binding variables
    
    # Set up dataframe to save info about selected eSNP-eGENE pairs 
    esnp_id <- paste0('eSNP_', id)
    es_id <- paste0('EffectSize_', id)
    snps_list <- snps$eSNP
    pairs <- data.frame(genes)
    pairs$eQTL <- NA
    pairs[, esnp_id] <- NA
    pairs[, es_id] <- 0
    genes2 <- data.frame(genes)
    
    for(i in 1:eqtl.n){
        again <- TRUE
        while (again == TRUE){
            if(length(snps_list) == 0) {
                stop("Not enough SNPs within desired MAF range. Increase the
                    eqtl.mafd allowed, include more SNPs, or reduce eqtl.n.")
            }
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
        
        pairs[pairs$geneID == match, ][, esnp_id] <- s
        pairs[pairs$geneID == match, ][, es_id] <- ES
        pairs[pairs$geneID == match, ]$eQTL <- 'global'
    }
    
    return(pairs)
}


#' Simulate group-specific eQTL associations. 
#' 
#' If the number of groups is >1, this function randomly selects eSNP-eGene
#' pairs to be group specific. The number of group specific pairs is controlled
#' by the eqtl.group.specific parameter. 
#'
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' @param pairs eGene:eSNP key dataframe.
#' @param groups array of group names
#' 
#' @return A dataframe eSNPs-eGenes pair assignments and their effect sizes 
#'
groupSpecificEffects <- function(eQTLparams, pairs, groups){
    
    g.specific.perc <- getParam(eQTLparams, "eqtl.group.specific")
    esnps <- unique(pairs$eSNP_g1)
    esnps <- esnps[!is.na(esnps)]
    g.spec <- sample(esnps, size=(length(esnps) * g.specific.perc))
    
    # Randomly sample the group assignment for group-specific eQTL 
    which.group <- sample(groups, length(g.spec), replace = T)
    which.genes <- pairs$geneID[pairs$eSNP_g1 %in% g.spec]
    
    for(g in groups){
        which.genes.g <- which.genes[which.group==g]
        not.g <- groups[groups != g]
        colNA <- names(pairs)[grepl(paste(paste0("eSNP_", not.g), collapse = "|"), 
                                    names(pairs) )]
        col0 <- names(pairs)[grepl(paste(paste0("EffectSize_", not.g), collapse = "|"),
                                   names(pairs) )]
        pairs[pairs$geneID %in% which.genes.g, colNA] <- NA
        pairs[pairs$geneID %in% which.genes.g, col0] <- 0.0
        pairs[pairs$geneID %in% which.genes.g, 'eQTL'] <- g
    }
    
    return(pairs)
}

#' Generate normalized mean gene expression matrix for whole eQTL population
#'
#' Use the approach outlined in Huang et. al 2018 (NAR) to assign normalized
#' mean expression levels for each gene for each sample. Where:
#' y = Effect Size * genotype + error, where error ~ Norm(0,1)
#'
#' @param id The group ID (e.g. "g1")
#' @param snps The dataframe with the genetic marker info
#' @param pairs A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' @param last_norm_matrix Dataframe of normalized means if they were already
#' simulated for another group.
#' @return normGeneMeansPop: matrix of normalized mean gene exp. levels.
#'
#' @importFrom stats rnorm
eQTLnormMeansMatrix <- function(id, snps, pairs, last_norm_matrix) {
    
    # Generate matrix of normalized mean expression values
    esnp_id <- paste0('eSNP_', id)
    es_id <- paste0('EffectSize_', id)
    samples <- names(snps)[grepl('Sample', names(snps))]
    
    # Make empty matrix or start from previous results
    if(id == 'g1'){
        genes <- pairs$geneID
        norm_matrix <- data.frame(matrix(ncol=length(samples),
                                         nrow=dim(pairs)[1], 
                                         dimnames=list(pairs$geneID, samples)))
    }else{
        norm_matrix <- last_norm_matrix
        genes <- subset(pairs, eQTL == id)$geneID
    }
    
    for(g in genes) { 
        
        ES <- pairs[pairs$geneID==g, es_id]
        error <-  rnorm(length(samples), mean = 0, sd = 1.5)
        
        if(ES == 0){
            norm_matrix[g,] <- error
            
        }else{
            eSNPsample <- pairs[pairs$geneID==g, esnp_id]
            genotype <- as.numeric(snps[snps$eSNP==eSNPsample, samples])
            norm_matrix[g,] <- (ES * genotype) + error
        }
    }
    
    return(norm_matrix)
}



#' Project normalized gene expression matrix into mean gene expression matrix 
#'
#' @param id The group ID (e.g. "g1")
#' @param pairs A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' @param nGeneMeansPop The normalized gene expression means for the population
#' @param last_MeansMatrix SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details.
#' 
#' @details
#' For each gene/sample, the normalized expression value (from rnorm) is 
#' transformed to the cumulative density function (pnorm) between 0 and 1, this
#' value is then quantile normalized (qnorm) using the gamma distribution 
#' parameterized from splatEstimate(). 
#'
#' @return MeansMatrix: Matrix of simulated gene means for eQTL population.
#' 
#' @importFrom stats pnorm qnorm sd quantile
#' @importFrom stats sd
#' 
eQTLMeansMatrix <- function(id, pairs, nGeneMeansPop, last_MeansMatrix){
    
    # Make empty matrix or start from previous results
    if(id == 'g1'){
        genes <- pairs$geneID
        MeansMatrix <- data.frame(nGeneMeansPop)
    }else{
        MeansMatrix <- last_MeansMatrix
        genes <- subset(pairs, pairs$eQTL == id)$geneID
    }
    
    # For each gene, simulate normal dist. of mean expression across population.
    for(g in genes){
        mean.gene <- pairs[pairs$geneID == g,]$expMean
        cv.gene <- pairs[pairs$geneID == g,]$expCV
        sd.gene <- cv.gene * mean.gene
        norm.mean <- mean(unlist(MeansMatrix[g, ]))
        norm.sd <- sd(unlist(MeansMatrix[g, ]))
        
        n_val <- as.numeric(MeansMatrix[g, ])
        pnorm.tmp <- pnorm(n_val, norm.mean, norm.sd)
        MeansMatrix[g, ] <- qnorm(pnorm.tmp, mean.gene, sd.gene)
        
    }
    MeansMatrix[MeansMatrix < 0] <- 0
    
    return(MeansMatrix)
}

    
#' Project normalized gene expression matrix into mean gene expression matrix 
#'
#' @param params SplatParams object containing parameters for the simulation.
#'               See \code{\link{SplatParams}} for details.
#' @param MeansMatrix The gene expression means for the population
#' 
#' @details
#' For each sample, expression value is quantile normalized (qgamma) using the
#' gamma distribution parameterized from splatEstimate(). 
#'
#' @return MeansMatrix: Matrix of simulated gene means for eQTL population.
#' 
#' @importFrom stats pnorm qgamma sd quantile

quantileNormalizeSC <- function(params, MeansMatrix){    
    
    # For each sample, quantile gamma normalize expression across genes.
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")

    for(s in names(MeansMatrix)){
        sample_mean <- mean(MeansMatrix[, s])
        sample_sd <- sd(MeansMatrix[, s])
        pnorm.tmp <- pnorm(MeansMatrix[, s], sample_mean, sample_sd)
        pnorm.tmp[pnorm.tmp == 1] <- 0.9999
        MeansMatrix[, s] <- qgamma(pnorm.tmp, shape=mean.shape, rate=mean.rate)
    }
    return(MeansMatrix)
}

