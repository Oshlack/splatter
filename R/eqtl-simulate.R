#' eQTL-Simulate
#'
#' Simulate mean gene counts for a population of samples based on sample
#' genotype with eQTL effects included for certain genes (eGenes).
#'
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details. Default=`newSplatParams()`.
#' @param eQTLparams eQTLParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details. Default=`neweQTLParams()`.
#' @param gff Data.frame of genes to simulate from a GFF or GTF file.
#' @param vcf Data.frame of genotypes for samples to simulate from a VCF file.
#' @param key Data.frame of complete or partial eQTL key that is output from
#'        `eQTLSimulate()`. If FALSE, a key will be generated from scratch.
#'        Default=FALSE.
#' @param verbose logical. Whether to print progress messages. Default=TRUE.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details SplatParams and eQTLParams can be set in a variety of ways. 1. If
#' not provided, default parameters are used. 2. Default parameters can be
#' overridden by supplying desired parameters using \code{\link{setParams}}.
#' 3. Parameters can be estimated from real data of your choice using
#' \code{\link{splatEstimate}} and \code{\link{eQTLEstimate}}.
#'
#' `eQTLSimulate` involves the following steps:
#' \enumerate{
#'     \item Load eQTL key or generate eQTL key from the GTF/GFF file.
#'     \item Format and subset genotype data from the VCF file.
#'     \item If not in key, assign expression mean and variance to each gene.
#'     \item If not in key, assign eGenes-eSNPs pairs and effect sizes.
#'     \item If not in key and groups >1, assign subset of eQTL associations as
#'     group-specific.
#'     \item Simulate mean gene expression matrix without eQTL effects
#'     \item Add eQTL effects to means matrix.
#'     \item Quantile normalize by sample to fit single-cell expression
#'     distribution as defined in `splatEstimate`.
#'     \item Add quantile normalized gene mean and cv info the eQTL key.}
#'
#' @return A list containing: `means` a dataframe or list of dataframes (if
#' n.groups > 1) with the simulated mean gene expression value for each gene
#' (row) and each sample (column) and `key` a dataframe with the eSNP-eGene
#' assignments.
#'
#' @seealso
#' \code{\link{eqtl.parse.gff}}, \code{\link{eqtl.parse.vcf}},
#' \code{\link{eqtl.assign.means}}, \code{\link{eqtl.assign.eqtl.effects}},
#' \code{\link{eqtl.assign.group.effects}}, \code{\link{eqtl.sim.means}},
#' \code{\link{eqtl.sim.eqtl.eff}}, \code{\link{eqtl.quan.norm.sc}},
#' \code{\link{eqtl.key.update.quan.norm}}

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

    if (verbose) {message("Loading VCF...")}
    vcf.parsed <- eqtl.parse.vcf(vcf, eQTLparams)
    snps <- vcf.parsed$snps
    samples <- vcf.parsed$samples
    groups <- paste0("g", seq(1, getParam(eQTLparams, "eqtl.groups")))

    # Read in genes and gene locations from GFF/GTF or from the provided key
    if (key == FALSE){
        if (verbose) {message("Pulling genes from GFF...")}
        key <- eqtl.parse.gff(gff)
    }else{
        if (verbose) {message("Using genes from key provided...")}
    }

    # If mean and CV not provided in key, simulate using params from eQTLparams
    if (!all(c("exp_mean", "exp_cv") %in% names(key))){
        if (verbose) {message("Assigning gene means & cv...")}
        key <- eqtl.assign.means(params, key, eQTLparams)
    }

    # If eqtl effects are not provided, simulate them.
    if (!all(c("eQTL", "eSNP", "EffectSize") %in% names(key))){
        if (verbose) {message("Assigning eQTL effects...")}
        key <- eqtl.assign.eqtl.effects(key, snps, eQTLparams)

        if(length(groups) > 1){
            key <- eqtl.assign.group.effects(key, groups, eQTLparams, params)}
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
        key <- eqtl.key.update.quan.norm(key, eMeansPopq_groups)

        if (verbose) {message("Done!")}
        return(list(means=eMeansPopq_groups, key=key))

    # Otherwise quantile normalize single matrix and return
    }else{
        eMeansPopq <- eqtl.quan.norm.sc(params, eMeansPop)
        key <- eqtl.key.update.quan.norm(key, eMeansPopq)

        if (verbose) {message("Done!")}
        return(list(means=eMeansPopq, key=key))
    }
}


#' Generate eQTL key matrix using information from the GTF/GFF file.
#'
#' @param gff Dataframe of genes to simulate from a GFF or GTF file.
#'
#' @return The partial eQTL key dataframe.
#'
eqtl.parse.gff <- function(gff){

    # Test input gff file
    if ((length(names(gff)) < 8 | nrow(gff[gff[,3]=="gene",]) < 1)) {
        stop("GFF file did not contain gene features or other issue with
            file format. See example data.")
    }

    genes <- gff[gff[,3]=="gene",]
    genes$geneID <- c(paste0("gene", 1:nrow(genes)))
    genes$chr <- genes[, 1]
    genes$gene_start <- genes[, 4]
    genes$gene_end <- genes[, 5]
    genes$gene_dir <- genes[, 7]
    genes$loc <- ifelse(genes$gene_dir == "-", genes$gene_end, genes$gene_start)
    key <- genes[,c("geneID", "chr", "gene_start", "gene_end", "gene_dir", "loc")]

    return(key)
}


#' Format and subset genotype data from a VCF file.
#'
#' Convert [0/0, 0/1, 1/1] alleles into [0,1,2] format and remove SNPs that do
#' not meet Minor Allele Frequency requirements specified in `eQTLParams`.
#'
#' @param vcf Data.frame of genotypes for samples to simulate from a VCF file.
#' @param eQTLparams eQTLparams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details. Default=`neweQTLParams()`.
#'
#' @return A list containing the genotype dataframe and a list of samples.
#'
eqtl.parse.vcf <- function(vcf, eQTLparams){

    eqtl.maf.min <- getParam(eQTLparams, "eqtl.maf.min")
    eqtl.maf.max <- getParam(eQTLparams, "eqtl.maf.max")

    # Read in genotype matrix in .vcf format
    vcf[, 3:9] <- NULL
    samp_n <- ncol(vcf) - 2
    samp_n_dig <- nchar(samp_n)
    if(names(vcf)[3] == "V10"){
        names(vcf) <- c("chr", "loc", paste0("Sample",
                                             formatC(1:samp_n, width=samp_n_dig,
                                                     format="d", flag="0")))
    }else{
        names(vcf)[1:2] <- c("chr", "loc")
    }

    samples <- names(vcf)[3:ncol(vcf)]
    vcf[] <- lapply(vcf, function(x) gsub("0/0", 0.0, x))
    vcf[] <- lapply(vcf, function(x) gsub("0/1", 1.0, x))
    vcf[] <- lapply(vcf, function(x) gsub("1/1", 2.0, x))
    vcf <- as.data.frame(sapply(vcf, as.numeric))
    vcf <- cbind(eSNP = paste("snp", vcf$chr, vcf$loc, sep=":"), vcf, stringsAsFactors=FALSE)

    # Filter out SNPs not within MAF requested
    vcf$MAF <- rowSums(vcf[, samples] / (length(samples) * 2))
    snps <- subset(vcf, MAF >= eqtl.maf.min & MAF <= eqtl.maf.max)
    row.names(snps) <- snps$eSNP
    snps$eSNP <- NULL

    return(list(snps=snps, samples=samples))
}


#' Assign expression mean and variance to each gene
#'
#' A mean and coefficient of variation is assigned to each gene by sampling from
#' gamma distributions parameterized from real data using `splatEstimate`.
#' The cv gamma distributions are binned by gene mean because the distribution
#' of variance in real data is not independent from the mean.
#'
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details.
#' @param key Partial eQTL key dataframe.
#' @param eQTLparams eQTLParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#'
#' @return The key updated with assigned means and variances.
#'
eqtl.assign.means <- function(params, key, eQTLparams){

    # Load parameters generated from real data using splatEstimate()
    pop_mean_shape <- getParam(eQTLparams, "pop.mean.shape")
    pop_mean_rate <- getParam(eQTLparams, "pop.mean.rate")
    cv.param <- getParam(eQTLparams, "pop.cv.param")

    # Sample gene means
    key$exp_mean <- rgamma(nrow(key), shape = pop_mean_shape,
                            rate = pop_mean_rate)
    key$exp_cv <- NULL

    # Sample coefficient of variation for each gene
    for (g in 1:nrow(key)){
        exp_mean <- key[g, "exp_mean"]
        bin <- cv.param[(cv.param$start < exp_mean) &
                            (cv.param$end >= exp_mean), ]
        key[g,"exp_cv"] <- rgamma(1, shape = bin$shape, rate = bin$rate)
    }

    return(key)
}


#' Assign eGenes-eSNPs pairs and effect sizes.
#'
#' Randomly pairs N genes (eGene) a SNP (eSNP) within the window size
#' (eqtl.dist) and assigns each pair an effect size sampled from a gamma
#' distribution parameterized using the effect sizes from a real eQTL study.
#'
#' @param key Partial eQTL key dataframe.
#' @param snps Dataframe of genotype information output from `vcf.parsed`.
#' @param eQTLparams eQTLParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details. Default=`neweQTLParams()`.
#'
#' @return The key updated with assigned eQTL effects.
#'
eqtl.assign.eqtl.effects <- function(key, snps, eQTLparams){
    eqtl.n <- getParam(eQTLparams, "eqtl.n")
    if (eqtl.n > nrow(key)){eqtl.n <- nrow(key)}
    eqtl.dist <- getParam(eQTLparams, "eqtl.dist")
    eqtlES_shape <- getParam(eQTLparams, "eqtl.ES.shape")
    eqtlES_rate <- getParam(eQTLparams, "eqtl.ES.rate")

    # Set up dataframe to save info about selected eSNP-eGENE pairs
    snps_list <- row.names(snps)
    key_tmp <- key
    key$eQTL <- NA
    key$eSNP <- NA
    key$EffectSize <- 0

    for(i in seq_len(eqtl.n)){
        again <- TRUE
        while (again == TRUE){
            if(length(snps_list) == 0) {
                stop("Not enough SNPs in MAF range.")
            }
            s <- sample(snps_list, 1)
            snps_list <- snps_list[!snps_list==s]

            l <- snps[s, "loc"]
            s_chr <- snps[s, "chr"]
            matches <- subset(key_tmp, (chr == s_chr & loc > l - eqtl.dist &
                                          loc < l + eqtl.dist))
            if(nrow(matches) > 0){
                match <- sample(matches$geneID, 1)
                again <- FALSE
            }
        }

        key_tmp <- key_tmp[!(key_tmp$geneID==match),]
        ES <- rgamma(1, shape = eqtlES_shape, rate = eqtlES_rate)

        key[key$geneID == match, ]$eSNP <- s
        key[key$geneID == match, ]$EffectSize <- ES
        key[key$geneID == match, ]$eQTL <- "global"

        # Randomly make some effects negative
        key$EffectSize <- key$EffectSize * sample(c(1, -1),
                                                  length(key$EffectSize),
                                                  replace = TRUE)
    }

    return(key)
}


#' Designate subset of eQTL associations as group-specific.
#'
#' If groups > 1, n eSNP-eGene pairs (n = 'eqtl.group.specific') are randomly
#' assigned as group specific.
#'
#' @param key Partial eQTL key dataframe.
#' @param groups array of group names
#' @param eQTLparams eQTLParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details. Default=`neweQTLParams()`.
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details. Default=`newSplatParams()`.
#'
#' @return he key updated with group eQTL and non-eQTL effects.
#'
eqtl.assign.group.effects <- function(key, groups, eQTLparams, params){

    # Assign group-specific eQTL
    eqtl.n <- getParam(eQTLparams, "eqtl.n")
    g.specific.perc <- getParam(eQTLparams, "eqtl.group.specific")
    n.groups <- length(groups)
    n.specific.each <- ceiling(eqtl.n * g.specific.perc / n.groups)

    for(g in groups){
        glob_genes <- subset(key, eQTL == "global")$geneID
        g.specific <- sample(glob_genes, size = n.specific.each)
        key$eQTL[key$geneID %in% g.specific] <- g
    }

    # Assign group-specific effects (differential expression, not eQTL)
    nGenes <- nrow(key)
    de.prob <- getParam(params, "de.prob")
    de.downProb <- getParam(params, "de.downProb")
    de.facLoc <- getParam(params, "de.facLoc")
    de.facScale <- getParam(params, "de.facScale")

    for (idx in seq_len(n.groups)) {
        de.facs <- getLNormFactors(nGenes, de.prob, de.downProb,
                                   de.facLoc, de.facScale)
        key[, paste0(groups[idx], "_group_effect")] <- de.facs
    }
    return(key)
}


#' Simulate mean gene expression matrix without eQTL effects
#'
#' Gene mean expression levels are assigned to each gene for each pair randomly
#' from a normal distribution parameterized using the mean and cv assigned to
#' each gene in the key.
#'
#' @param samples Vector containing the sample names.
#' @param key Partial eQTL key dataframe.
#'
#' @return Dataframe of gene mean expression levels WITHOUT eQTL effects.
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
#' Add eQTL effects and non-eQTL group effects to simulated means matrix.
#' The eQTL effects are incorporated using the following equation:
#' \deqn{Ygs = (ESg * Mg) + Mgs }
#' Where Ygs is the mean for gene *g* and sample *s*, ESg is the effect size
#' assigned to *g*, Mg is the mean expression assigned to *g*, and Mgs is the
#' mean expression sampled for *g* for *s*. Non-eQTL group effects are
#' incorporated as:
#' \deqn{Ygs = Mgs * GEg}
#' Where GEg is the group effect (i.e. differential expression) assigned to *g*.
#' To simulate multiple gene mean matrices with different group effects, this
#' function can be run with `id` designating the group id.
#'
#' @param id The group ID (e.g. "global" or "g1")
#' @param key Partial eQTL key dataframe.
#' @param snps Dataframe of genotype information output from `vcf.parsed`.
#' @param MeansPop Mean gene expression dataframe.
#'
#' @return Dataframe of gene mean expression levels WITH eQTL effects.
#'
eqtl.sim.eqtl.eff <- function(id, key, snps, MeansPop){

    # Add group-specific eQTL effects
    genes_use <- subset(key, eQTL == id)$geneID
    samples <- names(MeansPop)
    key$EffectSize_m <- key$exp_mean * key$EffectSize

    for(g in genes_use){
        without_eqtl <- as.numeric(MeansPop[g,])
        ES <- key[key$geneID == g, "EffectSize_m"]
        eSNPsample <- key[key$geneID == g, "eSNP"]
        genotype <- as.numeric(snps[eSNPsample, samples])
        MeansPop[g,] <- (ES * genotype) + without_eqtl
    }

    # Add group-specific non-eQTL effects
    if(id != "global"){
        MeansPop <- MeansPop * key[, paste0(id, "_group_effect")]
    }
    MeansPop[MeansPop < 0] <- 0

    return(MeansPop)
}


#' Quantile normalize by sample to fit sc expression distribution.
#'
#' For each sample, expression values are quantile normalized (qgamma)
#' using the gamma distribution parameterized from splatEstimate(). This ensures
#' the simulated gene means reflect the distribution expected from a sc dataset
#' and not a bulk dataset.
#'
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details. Default=`newSplatParams()`.
#' @param MeansMatrix Mean gene expression with eQTL effects dataframe.
#'
#' @return Dataframe of quantile normalized gene mean expression levels.
#'
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @export

eqtl.quan.norm.sc <- function(params, MeansMatrix){

    # Generate sample target distribution from sc parameters
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    target <- rgamma(10000, shape=mean.shape, rate=mean.rate)

    mat_norm <- preprocessCore::normalize.quantiles.use.target(as.matrix(MeansMatrix), target)
    mat_norm[mat_norm < 0] <- 0
    df_norm <- as.data.frame(mat_norm, row.names = row.names(MeansMatrix))
    names(df_norm) <- names(MeansMatrix)

    return(df_norm)
}

#' Add quantile normalized gene mean and cv info the eQTL key.
#'
#' @param key Partial eQTL key dataframe.
#' @param MeansMatrix The output from `eqtl.quan.norm.sc` as a matrix or list of
#'                    matrices.
#'
#' @return Final eQTL key.
#'
#' @export
eqtl.key.update.quan.norm <- function(key, MeansMatrix){

    if (type(MeansMatrix) == "list"){
        qn_means <- list()
        qn_cvs <- list()

        for(group in names(MeansMatrix)){
            qn_mean_gr <- rowMeans(MeansMatrix[[group]])
            qn_cv_gr <- apply(MeansMatrix[[group]], 1, FUN=co.var)
            qn_means[[group]] <- qn_mean_gr
            qn_cvs[[group]] <- qn_cv_gr
        }
        qn_mean <- rowMeans(as.data.frame(qn_means))
        qn_cv <- rowMeans(as.data.frame(qn_cvs))

    }else{
        qn_mean <- rowMeans(MeansMatrix)
        qn_cv <- apply(MeansMatrix, 1, FUN=co.var)
    }

    qn_df <- data.frame(list(expQN_mean=qn_mean, expQN_cv=qn_cv))
    qn_df$geneID <- row.names(qn_df)
    key <- merge(key, qn_df, by="geneID")

    return(key)
}
