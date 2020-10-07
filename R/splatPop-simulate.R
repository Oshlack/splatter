#' splatPop simulation
#'
#' Simulate scRNA-seq count data using the splat model for a population of
#' individuals with correlation structure.
#'
#' @param params splatPopParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#'
#' This functions is for simulating data in a single step. It consists of a
#' call to \code{\link{splatPopMeans}}, which simulates a mean expression level
#' per gene per sample, followed by a call to \code{\link{splatPopSC}}, which
#' uses the splat model to simulate single-cell counts per individual. Please
#' see the documentation for those functions for more details.
#'
#' @seealso
#' \code{\link{splatPopMeans}}, \code{\link{splatPopSC}}
#'
#' @return SingleCellExperiment containing simulated counts and intermediate
#' values
#'
#' @examples
#'
#' sim <- splatPopSimulate()
#'
#' @export
splatPopSimulate <- function(params = newsplatPopParams(), verbose = TRUE,
                             ...) {
    
    sim_means <- splatPopMeans(params, verbose, ...)
    sim_sc <- splatPopSC(sim_means, params, verbose)
    
    return(sim)
}


#' popSimulate
#'
#' Simulate mean gene counts for a population of samples based on sample
#' genotype with eQTL effects included for certain genes (eGenes).
#'
#' @param params splatPopParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{splatPopParams}} for details. Default=`newSplatPopParams()`.
#' @param genes Specify genes to include. Either provide data.frame of genes
#'        from a GFF/GTF or set to "random" to have popSimulate randomly
#'        generate genes. Default="random".
#' @param vcf Data.frame of genotypes for samples to simulate from a VCF file.
#' @param key Data.frame of complete or partial eQTL key that is output from
#'        `popSimulate()`. If FALSE, a key will be generated from scratch.
#'        Default=FALSE.
#' @param verbose logical. Whether to print progress messages. Default=TRUE.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details splatPopParams can be set in a variety of ways. 1. If
#' not provided, default parameters are used. 2. Default parameters can be
#' overridden by supplying desired parameters using \code{\link{setParams}}.
#' 3. Parameters can be estimated from real data of your choice using
#' \code{\link{splatEstimate}} and \code{\link{popEstimate}}.
#'
#' `popSimulate` involves the following steps:
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
#' \code{\link{pop.parse.vcf}}, \code{\link{pop.parse.gff}},
#' \code{\link{pop.random.genes}}, \code{\link{pop.assign.means}},
#' \code{\link{pop.assign.eqtl.effects}},
#' \code{\link{pop.assign.group.effects}}, \code{\link{pop.sim.means}},
#' \code{\link{pop.sim.eqtl.eff}}, \code{\link{pop.quan.norm.sc}},
#' \code{\link{pop.key.update.quan.norm}}
#'
#' @export
#'

splatPopMeans <- function(vcf, params = newSplatPopParams(), verbose = TRUE, ...) {
    
    checkmate::assertClass(params, "splatPopParams")
    params <- setParams(params, ...)
    
    # Set random seed
    seed <- getParam(splatPopParams, "seed")
    set.seed(seed)
    
    # Read in SNPs from VCF
    if (verbose) {message("Loading VCF...")}
    vcf.parsed <- pop.parse.vcf(vcf, params)
    snps <- vcf.parsed$snps
    samples <- vcf.parsed$samples
    groups <- paste0("g", seq(1, getParam(params, "eqtl.groups")))
    
    # Read in genes and gene locations from GFF/GTF or from the provided key
    if (key == FALSE){
        if(genes == "random"){
            if (verbose) {message("Generating random genes...")}
            key <- pop.random.genes(params, vcf)
        }else{
            if (verbose) {message("Pulling genes from GFF...")}
            key <- pop.parse.gff(genes)
        }
        
    }else{
        if (verbose) {message("Using genes from key provided...")}
    }
    
    # If mean and CV not provided in key, simulate using params from params
    if (!all(c("exp_mean", "exp_cv") %in% names(key))){
        if (verbose) {message("Assigning gene means & cv...")}
        key <- pop.assign.means(params, key, params)
    }
    
    # If eqtl effects are not provided, simulate them.
    if (!all(c("eQTL", "eSNP", "EffectSize") %in% names(key))){
        if (verbose) {message("Assigning eQTL effects...")}
        key <- pop.assign.eqtl.effects(key, snps, params)
        
        if(length(groups) > 1){
            key <- pop.assign.group.effects(key, groups, params, params)}
    }
    
    # Use info in key to generate matrix of gene expression levels
    if (verbose) {message("Simulating gene means for population...")}
    MeansPop <- pop.sim.means(samples, key)
    MeansPop <- pop.quan.norm.sc(params, MeansPop) # keep if switching quant.norm order
    eMeansPop <- pop.sim.eqtl.eff('global', key, snps, MeansPop)
    
    # Add group-specific eQTL effects if present
    if(length(groups) > 1){
        eMeansPopq_groups <- list()
        
        for(id in groups){
            eMeansPop.g <- pop.sim.eqtl.eff(id, key, snps, eMeansPop)
            #eMeansPop.g.q <- pop.quan.norm.sc(params, eMeansPop.g)
            #eMeansPopq_groups[[id]] <- eMeansPop.g.q
            eMeansPop.g[eMeansPop.g <= 0] <- 1e-4 # keep if switching quant.norm order
            eMeansPopq_groups[[id]] <- eMeansPop.g # keep if switching quant.norm order
            
        }
        key <- pop.key.update.quan.norm(key, eMeansPopq_groups)
        
        if (verbose) {message("Done!")}
        return(list(means=eMeansPopq_groups, key=key))
        
        # Otherwise quantile normalize single matrix and return
    }else{
        #eMeansPopq <- pop.quan.norm.sc(params, eMeansPop)
        #key <- pop.key.update.quan.norm(key, eMeansPopq)
        eMeansPop[eMeansPop <= 0] <- 1e-4 # keep if switching quant.norm order
        key <- pop.key.update.quan.norm(key, eMeansPop) # keep if switching quant.norm order
        if (verbose) {message("Done!")}
        return(list(means=eMeansPop, key=key))
        #return(list(means=eMeansPopq, key=key)) # keep if switching quant.norm order
    }
}


#' Format and subset genotype data from a VCF file.
#'
#' Convert [0/0, 0/1, 1/1] alleles into [0,1,2] format and remove SNPs that do
#' not meet Minor Allele Frequency requirements specified in `splatPopParams`.
#'
#' @param vcf Data.frame of genotypes for samples to simulate from a VCF file.
#' @param params splatPopParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{splatPopParams}} for details. Default=`newSplatPopParams()`.
#'
#' @return A list containing the genotype dataframe and a list of samples.
#'
pop.parse.vcf <- function(vcf, params){
    
    eqtl.maf.min <- getParam(params, "eqtl.maf.min")
    eqtl.maf.max <- getParam(params, "eqtl.maf.max")
    
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
splatPop.parse.vcf <- function(vcf, params){
    
    eqtl.maf.min <- getParam(params, "eqtl.maf.min")
    eqtl.maf.max <- getParam(params, "eqtl.maf.max")
    
    # Read in genotype matrix in .vcf format
    vcf_object <- vcf@gt
    vcf_loc <- vcf@fix
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

#' Generating eQTL key matrix with random genes
#'
#' @param params splatPopParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{splatPopParams}} for details. Default=`newSplatPopParams()`.
#' @param vcf output from `pop.parse.vcf`
#'
#' @return The partial eQTL key dataframe.
#'
pop.random.genes <- function(params, vcf){
    
    nGenes <- getParam(params, "nGenes")
    
    gene_numbers <- sprintf(paste0("%0", nchar(nGenes), "d"), 1:nGenes)
    key <- data.frame(list('geneID' = paste0('gene', gene_numbers),
                           'tmp_snp' = sample(1:nrow(vcf), nGenes)))
    key$chr <- vcf[key$tmp_snp, 'V1']
    key$loc <- vcf[key$tmp_snp, 'V2'] + sample(-1e3:1e3, nGenes, replace=TRUE)
    
    key <- key[,c("geneID", "chr", "loc")]
    
    return(key)
}

#' Generate eQTL key matrix using information from the GTF/GFF file.
#'
#' @param gff Dataframe of genes to simulate from a GFF or GTF file.
#'
#' @return The partial eQTL key dataframe.
#'
pop.parse.gff <- function(gff){
    
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



#' Assign expression mean and variance to each gene
#'
#' A mean and coefficient of variation is assigned to each gene by sampling from
#' gamma distributions parameterized from real data using `splatEstimate`.
#' The cv gamma distributions are binned by gene mean because the distribution
#' of variance in real data is not independent from the mean.
#'
#' @param params splatPopParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{splatPopParams}} for details.
#' @param key Partial eQTL key dataframe.
#'
#' @return The key updated with assigned means and variances.
#'
pop.assign.means <- function(params, key){
    
    # Load parameters generated from real data using splatEstimate()
    pop_mean_shape <- getParam(params, "pop.mean.shape")
    pop_mean_rate <- getParam(params, "pop.mean.rate")
    cv.param <- getParam(params, "pop.cv.param")
    
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
#' @param params splatPopParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{splatPopParams}} for details. Default=`newSplatPopParams()`.
#'
#' @return The key updated with assigned eQTL effects.
#'
pop.assign.eqtl.effects <- function(key, snps, params){
    
    eqtl.n <- getParam(params, "eqtl.n")
    if (eqtl.n > nrow(key)){eqtl.n <- nrow(key)} # Can't be greater than nGenes
    if (eqtl.n <= 1){eqtl.n <- nrow(key) * eqtl.n} # If <= 1 it is a percent
    
    eqtl.dist <- getParam(params, "eqtl.dist")
    eqtlES_shape <- getParam(params, "eqtl.ES.shape")
    eqtlES_rate <- getParam(params, "eqtl.ES.rate")
    
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
#' @param params splatPopParams object containing parameters for the
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{splatPopParams}} for details. Default=`newSplatPopParams()`.
#'
#' @return he key updated with group eQTL and non-eQTL effects.
#'
pop.assign.group.effects <- function(key, groups, params){
    
    # Assign group-specific eQTL
    eqtl.n <- getParam(params, "eqtl.n")
    g.specific.perc <- getParam(params, "eqtl.group.specific")
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
pop.sim.means <- function(samples, key){
    
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
pop.sim.eqtl.eff <- function(id, key, snps, MeansPop){
    
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

pop.quan.norm.sc <- function(params, MeansMatrix){
    
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
#' @param MeansMatrix The output from `pop.quan.norm.sc` as a matrix or list of
#'                    matrices.
#'
#' @return Final eQTL key.
#'
#' @export
pop.key.update.quan.norm <- function(key, MeansMatrix){
    
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
