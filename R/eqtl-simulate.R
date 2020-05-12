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
#' @details Parameters can be set in a variety of ways. If no parameters are
#' provided the default parameters are used. Any parameters in \code{params} or
#' \code{eQTLparams} can be overridden by supplying additional arguments through
#' a call to \code{\link{setParams}}. This design allows the user flexibility in
#' how they supply parameters and allows small adjustments without creating a
#' new \code{SplatParams} object. See examples for a demonstration of how this
#' can be used.
#'
#' eQTLSimulate involves the following steps:
#' \enumerate{
#'     \item Read and format genotype data given a .vcf.
#'     \item Read in genes from GTF/GFF or from a provided key.
#'     \item If not provided in key, assign expression mean and cv to each gene.
#'     \item If not provided in key, assign eQTL and effect sizes to n genes.
#'     \item Assign gene expression values for each individual for each gene
#'     by sampling randomly from a normal distribution parameterized from
#'     eqtlEstimate.
#'     \item Add in eQTL effects.
#'     \item Quantile normalize expression levels by sample from a gamma
#'     distribution parameterized from splatEstimate.}
#'
#' @return A list continaing: `means` a dataframe or list of dataframes (if
#' n.groups > 1) with the simulated mean gene expression value for each gene
#' (row) and each sample (column) and `key` a dataframe with the eSNP-eGene
#' assignments.
#'
#' @seealso
#' \code{\link{eqtl.parse.vcf}}, \code{\link{eqtl.parse.snps}},
#' \code{\link{eqtl.assign.means}}, \code{\link{eqtl.assign.eqtl.effects}},
#' \code{\link{eqtl.assign.group.effects}}, \code{\link{eqtl.sim.means}}
#' \code{\link{eqtl.sim.eqtl.eff}}, \code{\link{eqtl.quan.norm.sc}}

#' @examples
#' # Load example data
#' data(ex_gff)
#' data(ex_snps)
#' pop.sim <- eQTLSimulate()
#'
#' @export
#' @importFrom utils read.csv
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

    if (verbose) {message("Loading VCF...")}
    snps <- eqtl.parse.snps(vcf, eQTLparams)
    samples <- names(snps)[grepl('Sample', names(snps))]
    groups <- paste0('g', seq(1, getParam(eQTLparams, "eqtl.groups")))

    # Read in genes and gene locations from GFF/GTF or from the provided key
    if (key == FALSE){
        key <- eqtl.parse.vcf(gff)
    }else{
        key <- read.csv(key)
    }

    # If mean and CV not provided in key, sim using parameters from eQTLparams.
    if (!all(c("exp_mean", "exp_cv") %in% names(key))){
        if (verbose) {message("Assigning gene means & cv...")}
        key <- eqtl.assign.means(params, key, eQTLparams)
    }

    # If eqtl effects are not provided, simulate them.
    if (!all(c("eQTL", "eSNP", "EffectSize") %in% names(key))){
        if (verbose) {message("Assigning eQTL effects...")}
        key <- eqtl.assign.eqtl.effects(key, snps, eQTLparams)
        if(length(groups) > 1){key <- eqtl.assign.group.effects(key, groups, eQTLparams)}
    }

    # Use info in key to generate matrix of gene expression levels
    if (verbose) {message("Simulating gene means for population...")}
    MeansPop <- eqtl.sim.means(samples, key)
    eMeansPop <- eqtl.sim.eqtl.eff('global', key, snps, MeansPop)

    # Add group-specific eQTL effects if present
    if(length(groups) > 1){
        eMeansPopq_groups <- list()
        for(id in groups){
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
#' Select features labeled "gene", pull the chromosome and transcriptional start
#' site (considering strand direction) for each gene.
#'
#' @param gff Dataframe containing the genes to include in GFF/GTF format.
#'
#' @return A dataframe containing gene IDs, chr, and location.
#'
eqtl.parse.vcf <- function(gff){

    # Test input gff file
    if ((length(names(gff)) < 8 | nrow(gff[gff[,3]=="gene",]) < 1)) {
        stop("GFF file did not contain gene features or other issue with
            file format. See example data.")
    }

    genes <- gff[gff[,3]=="gene",]
    genes$loc <- ifelse(genes[, 7] == "-", genes[, 5], genes[, 4])
    genes$geneID <- c(paste0("gene", 1:dim(genes)[1]))
    genes$chr <- genes[, 1]
    key <- genes[,c('geneID', 'chr', 'loc')]

    return(key)
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

    # Test input VCF file
    if (!('V1' %in% names(vcf))) {
        stop("snps not in the expected VCF format. See example data.")
    }

    # Read in genotype matrix in .vcf format
    vcf[, c('V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9')] <- NULL
    names(vcf) <- c('chr', 'loc', paste0('Sample', 1:(dim(vcf)[2]-2)))

    vcf[] <- lapply(vcf, function(x) gsub("0/0", 0.0, x))
    vcf[] <- lapply(vcf, function(x) gsub("0/1", 1.0, x))
    vcf[] <- lapply(vcf, function(x) gsub("1/1", 2.0, x))
    vcf <- as.data.frame(sapply(vcf, as.numeric))
    vcf <- cbind(eSNP = paste('snp', vcf$chr, vcf$loc, sep=":"), vcf, stringsAsFactors=FALSE)

    # Filter out SNPs not within MAF requested
    samples <- names(vcf)[grepl('Sample', names(vcf))]
    vcf$MAF <- rowSums(vcf[, samples] / (length(samples) * 2))
    snps <- subset(vcf, MAF > eqtl.maf-eqtl.mafd & MAF < eqtl.maf+eqtl.mafd)
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
eqtl.assign.means <- function(params, genes, eQTLparams){

    # Load parameters generated from real data using splatEstimate()
    bulk_mean_shape <- getParam(eQTLparams, "bulkmean.shape")
    bulk_mean_rate <- getParam(eQTLparams, "bulkmean.rate")
    cv.param <- getParam(eQTLparams, "bulkcv.param")

    # Sample gene means
    genes$exp_mean <- rgamma(nrow(genes), shape = bulk_mean_shape,
                            rate = bulk_mean_rate)
    genes$exp_cv <- NULL

    # Sample coefficient of variation for each gene
    for (g in 1:nrow(genes)){
        mean <- genes[g, 'exp_mean']
        bin <- cv.param[(cv.param$start < mean) & (cv.param$end >= mean), ]
        genes[g,'exp_cv'] <- rgamma(1, shape = bin$shape, rate = bin$rate)
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
#' @param genes Dataframe with gene ID, chromosome, chromosome, and location
#' @param snps Dataframe with SNP ID, location, and sample genotypes
#' @param eQTLparams eQTLParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#'
#' @return A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' @importFrom dplyr mutate "%>%"
#'
eqtl.assign.eqtl.effects <- function(genes, snps, eQTLparams){
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
            s_chr <- snps[s, 'chr']
            matches <- subset(genes, (chr == s_chr & loc > l - eqtl.dist &
                                          loc < l + eqtl.dist))
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
eqtl.assign.group.effects <- function(key, groups, eQTLparams){

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
#' each gene (see `eqtl.assign.means`).
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
                                      mean = key[key$geneID == g,]$exp_mean,
                                      sd = key[key$geneID == g,]$exp_mean *
                                          key[key$geneID == g,]$exp_cv))

    means.df <- data.frame(do.call(rbind, means), row.names = key$geneID)
    names(means.df) <- samples

    return(means.df)
}


#' Add eQTL effects to means matrix
#'
#' For eSNP-eGene pairs for id group, the eQTL effect is added to the
#' gene mean y = (scaled Effect Size) \* (assigned mean for that gene)
#' \* genotype + error
#'
#' @param id The group ID (e.g. "global" or "g1")
#' @param key A dataframe eSNPs-eGenes pair assignments and their effect sizes
#' @param snps A dataframe containing SNP names, locations, & sample genotypes.
#' @param MeansPop Mean gene expression dataframe
#'
#' @return Matrix of simulated gene means for eQTL population.
#'
eqtl.sim.eqtl.eff <- function(id, key, snps, MeansPop){

    # Get list of genes with eQTL effects of the type specified in id
    genes_use <- subset(key, eQTL == id)$geneID
    samples <- names(snps)[grepl('Sample', names(snps))]

    # Calculate eSNP effect size given gean mean assigned
    key$EffectSize_m <- key$exp_mean * key$EffectSize

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
#' For each sample, expression value is quantile normalized (qgamma)
#' using the gamma distribution parameterized from splatEstimate().
#'
#' @param params SplatParams object containing parameters for the simulation.
#'               See \code{\link{SplatParams}} for details.
#' @param MeansMatrix The gene expression means for the population
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
