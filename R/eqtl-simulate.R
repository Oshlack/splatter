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
#' @param key Path to file with eQTL key information. If FALSE (default), will 
#'        generate new key using params in eQTLparams and genes from gff
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
#'     \item (optional) Save eQTL key
#' }
#'
#' @return A list continaing: `means` a dataframe or list of dataframes (if 
#' n.groups > 1) with the simulated mean gene expression value for each gene 
#' (row) and each sample (column) and `key` a dataframe with the eSNP-eGene
#' assignments.
#' 
#' @seealso
#' \code{\link{eqtl.parse.vcf}}, \code{\link{eqtl.parse.snps}},
#' \code{\link{eqtl.sample.means}}, \code{\link{eqtl.key}}, 
#' \code{\link{eqtl.group.effects}}, \code{\link{eqtl.sim.means}}
#' \code{\link{eqtl.sim.eqtl.eff}}, \code{\link{eqtl.quan.norm.sc}}

#' @examples
#' # Load example data
#' data(ex_gff)
#' data(ex_snps)
#' pop.sim <- eQTLSimulate()
#' 
#' @export
#' 
eQTLSimulate <- function(params = newSplatParams(),
                      eQTLparams = neweQTLParams(),
                      gff = ex_gff, 
                      vcf = ex_snps,
                      key = FALSE,
                      verbose = TRUE, ...) {
    
    # Set random seed
    seed <- getParam(eQTLparams, "seed")
    set.seed(seed)

    if (verbose) {message("Loading SNP data...")}
    snps <- eqtl.parse.snps(vcf, eQTLparams)
    samples <- names(snps)[grepl('Sample', names(snps))]
    groups <- paste0('g', seq(1, getParam(eQTLparams, "eqtl.groups")))
    
    if (key == FALSE){
        if (verbose) {message("Generating eGene-eSNP key...")}
        genes <- eqtl.parse.vcf(gff)
        genes <- eqtl.sample.means(params, genes, eQTLparams)
        key <- eqtl.key(genes, snps, eQTLparams)
        if(length(groups) > 1){key <- eqtl.group.effects(key, groups, eQTLparams)}
        
    }else{
        key <- read.csv(key)
    }
    
    if (verbose) {message("Simulating gene means with global eQTL effects...")}
    MeansPop <- eqtl.sim.means(samples, key)
    eMeansPop <- eqtl.sim.eqtl.eff('global', key, snps, MeansPop)
    
    # Add group-specific eQTL effects if present
    if(length(groups) > 1){
        
        eMeansPopq_groups <- list()
        for(id in groups){
            
            if (verbose) {message(paste0("Adding ", id, "-specific eQTL effects..."))}
            eMeansPop.g <- eqtl.sim.eqtl.eff(id, key, snps, eMeansPop)
            eMeansPop.g.q <- eqtl.quan.norm.sc(params, eMeansPop.g)
            eMeansPopq_groups[[id]] <- eMeansPop.g.q
        }
        
        if (verbose) {message("Done!")} 
        return(list(means=eMeansPopq_groups, key=key))
    
    # Otherwise quantile normalize single matrix and return
    }else{
        
        eMeansPopq <- eqtl.quan.norm.sc(params, eMeansPop)
        
        if (verbose) {message("Done!")} 
        return(list(means=eMeansPopq, key=key))
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
#' 
eqtl.parse.vcf <- function(gff){
    
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
#' 
eqtl.parse.snps <- function(vcf, eQTLparams){
    
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
    row.names(snps) <- snps$eSNP
    snps$eSNP <- NULL
    
    return(snps)
}


#' Assign a mean expression value to each gene
#'
#' Means are sampled from a gamma distribution parameterized by splatEstimate, 
#' then sample a coefficient of variation for each gene from a gamma 
#' distribution parameterized by splatEstimate for genes in expression bins.
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
#' 
eqtl.sample.means <- function(params, genes, eQTLparams){
    
    # Load parameters generated from real data using splatEstimate()
    bulk_mean_shape <- getParam(eQTLparams, "bulkmean.shape")
    bulk_mean_rate <- getParam(eQTLparams, "bulkmean.rate")
    cv.param <- getParam(eQTLparams, "bulkcv.param")
    
    # Sample gene means
    genes$expMean <- rgamma(nrow(genes), shape = bulk_mean_shape, 
                            rate = bulk_mean_rate)
    genes$expCV <- NULL
    
    # Sample coefficient of variation for each gene
    for (g in 1:nrow(genes)){
        mean <- genes[g, 'expMean']
        bin <- cv.param[(cv.param$start < mean) & (cv.param$end >= mean), ]
        genes[g,'expCV'] <- rgamma(1, shape = bin$shape, rate = bin$rate)
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
#' @param genes Dataframe with gene ID and location
#' @param snps Dataframe with SNP ID, location, and sample genotypes
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' 
#' @return A dataframe eSNPs-eGenes pair assignments and their effect sizes 
#' @importFrom dplyr mutate "%>%"
#' 
eqtl.key <- function(genes, snps, eQTLparams){
    eqtl.n <- getParam(eQTLparams, "eqtl.n")
    if (eqtl.n > dim(genes)[1]){
        eqtl.n <- dim(genes)[1]
    }
    eqtl.dist <- getParam(eQTLparams, "eqtl.dist")
    eqtlES_shape <- getParam(eQTLparams, "eqtlES.shape") 
    eqtlES_rate <- getParam(eQTLparams, "eqtlES.rate")

    # Set up dataframe to save info about selected eSNP-eGENE pairs 
    snps_list <- row.names(snps)
    key <- genes %>% mutate(eQTL = NA, eSNP = NA, EffectSize = 0)

    for(i in 1:eqtl.n){
        again <- TRUE
        while (again == TRUE){
            if(length(snps_list) == 0) {
                stop("Not enough SNPs within desired MAF range. Increase the
                    eqtl.mafd allowed, include more SNPs, or reduce eqtl.n.")
            }
            s <- sample(snps_list, 1)
            snps_list <- snps_list[!snps_list==s]
            
            l <- snps[s, 'loc']
            matches <- subset(genes, TSS > l - eqtl.dist & TSS < l + eqtl.dist)
            if(dim(matches)[1] > 0){
                match <- sample(matches$geneID, 1)
                again <- FALSE
            }
        }
        
        genes <- genes[!(genes$geneID==match),]
        ES <- rgamma(1, shape = eqtlES_shape, rate = eqtlES_rate)
        
        key[key$geneID == match, ]$eSNP <- s
        key[key$geneID == match, ]$EffectSize <- ES
        key[key$geneID == match, ]$eQTL <- 'global'
        
        # Randomly make some effects negative
        key$EffectSize <- key$EffectSize * sample(c(1, -1), 
                                                  length(key$EffectSize), 
                                                  replace = TRUE)
    }
    
    return(key)
}


#' Simulate group-specific eQTL associations. 
#' 
#' If the number of groups is >1, this function randomly selects eSNP-eGene
#' pairs to be group specific. The number of group specific pairs is controlled
#' by the eqtl.group.specific parameter. 
#'
#' @param key eGene:eSNP key dataframe.
#' @param groups array of group names
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' 
#' @return A dataframe eSNPs-eGenes pair assignments and their effect sizes 
#'
eqtl.group.effects <- function(key, groups, eQTLparams){
    
    eqtl.n <- getParam(eQTLparams, "eqtl.n")
    g.specific.perc <- getParam(eQTLparams, "eqtl.group.specific")
    n.groups <- length(groups)
    n.specific.each <- ceiling(eqtl.n * g.specific.perc / n.groups)
    
    for(g in groups){
        glob_genes <- subset(key, eQTL == 'global')$geneID
        g.specific <- sample(glob_genes, size = n.specific.each)
        key$eQTL[key$geneID %in% g.specific] <- g
    }
    
    return(key)
}


#' Simulate mean gene expression matrix without eQTL effects
#'
#' Gene mean expression levels are assigned to each gene:sample pair randomly
#' from a normal distribution parameterized using the mean and CV assigned to 
#' each gene (see `eqtl.sample.means`). 
#'
#' @param samples Vector with the names of the samples
#' @param key A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' 
#' @return Mean gene expression levels for the population WITHOUT eQTL effects
#'
#' @importFrom stats rnorm
#' 
eqtl.sim.means <- function(samples, key){
    
    means <- lapply(key$geneID, 
                    function(g) rnorm(length(samples), 
                                      mean = key[key$geneID == g,]$expMean,
                                      sd = key[key$geneID == g,]$expMean *
                                          key[key$geneID == g,]$expCV))
    
    means.df <- data.frame(do.call(rbind, means), row.names = key$geneID)
    names(means.df) <- samples

    return(means.df)
}
    

#' Add eQTL effects to means matrix
#'
#' @param id The group ID (e.g. "global" or "g1")
#' @param key A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' @param snps A dataframe containing SNP names, locations, & sample genotypes.
#' @param MeansPop Mean gene expression dataframe
#' 
#' @details
#' For eSNP-eGene pairs for id group, the eQTL effect is added to the gene mean
#' y = (scaled Effect Size) \* (assigned mean for that gene) \* genotype + error
#'
#' @return Matrix of simulated gene means for eQTL population.
#' 
eqtl.sim.eqtl.eff <- function(id, key, snps, MeansPop){

    # Get list of genes with eQTL effects of the type specified in id
    genes_use <- subset(key, eQTL == id)$geneID
    samples <- names(snps)[grepl('Sample', names(snps))]
    
    # Calculate eSNP effect size given gean mean assigned
    key$EffectSize_m <- key$expMean * key$EffectSize
    
    for(g in genes_use){
        without_eqtl <- as.numeric(MeansPop[g,])
        ES <- subset(key, key$geneID == g)$EffectSize_m
        eSNPsample <- subset(key, key$geneID == g)$eSNP
        genotype <- as.numeric(snps[eSNPsample, samples])
        MeansPop[g,] <- (ES * genotype) + without_eqtl
    }
    
    MeansPop[MeansPop < 0] <- 0
    
    return(MeansPop)
}

#' Quantile normalize expression levels by sample to fit sc parameters  
#'
#' @param params SplatParams object containing parameters for the simulation.
#'               See \code{\link{SplatParams}} for details.
#' @param MeansMatrix The gene expression means for the population
#' 
#' @details
#' For each sample, expression value is quantile normalized (qgamma) using the
#' gamma distribution parameterized from splatEstimate(). 
#'
#' @return Matrix of simulated gene means for eQTL population.
#' 
#' @importFrom stats pnorm qgamma sd quantile

eqtl.quan.norm.sc <- function(params, MeansMatrix){    
    
    # For each sample, quantile gamma normalize expression across genes.
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")

    for(s in names(MeansMatrix)){
        s_values <- MeansMatrix[, s]
        pnorm.tmp <- pnorm(s_values, mean(s_values), sd(s_values))
        pnorm.tmp[pnorm.tmp == 1] <- 1 - 1e-6
        MeansMatrix[, s] <- qgamma(pnorm.tmp, shape=mean.shape, rate=mean.rate)
    }
    
    MeansMatrix[MeansMatrix < 0] <- 0
    
    return(MeansMatrix)
}

