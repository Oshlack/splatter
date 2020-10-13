#' splatPop simulation
#'
#' Simulate scRNA-seq count data using the splat model for a population of
#' individuals with correlation structure.
#'
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details.
#' @param vcf vcfR object containing genotypes of samples to simulate.
#' @param method which simulation method to use. Options are "single" which
#'        produces a single population, "groups" which produces distinct groups
#'        (eg. cell types), "paths" which selects cells from continuous
#'        trajectories (eg. differentiation processes).
#' @param gff Either NULL or a data.frame object containing a GFF/GTF file.
#' @param key Either NULL or a data.frame object containing a full or partial
#'        splatPop key.
#' @param counts.only logical. Whether to save only counts in sce object.       
#' @param verbose logical. Whether to print progress messages. 
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#'
#' This functions is for simulating data in a single step. It consists of a
#' call to \code{\link{splatPopSimulateMeans}}, which simulates a mean 
#' expression level per gene per sample, followed by a call to 
#' \code{\link{splatPopSimulateSC}}, which uses the splat model to simulate 
#' single-cell counts per individual. Please see the documentation for those 
#' functions for more details.
#'
#' @seealso
#' \code{\link{splatPopSimulateMeans}}, \code{\link{splatPopSimulateSC}}
#'
#' @return A list containing: `means` a data.frame (or list of data.frames if
#' n.groups > 1) with the simulated mean gene expression value for each gene
#' (row) and each sample (column), `key` a data.frame with population information
#' including eQTL and group effects, and `sce` a SingleCellExperiment object
#' containing simulated counts and intermediate values.
#'
#' @examples
#' if (requireNamespace(c("vcfR", "preprocessCore"), quietly = TRUE)) {
#'     library(vcfR)
#'     vcf <- mockVCF(n.samples = 6)
#'     sim <- splatPopSimulate(vcf = vcf)
#' }
#' sim <- splatPopSimulate()
#'
#' @export
splatPopSimulate <- function(params = newSplatPopParams(nGenes = 1000),
                             vcf = mockVCF(), 
                             method = c("single", "groups", "paths"),
                             gff = NULL,
                             key = NULL,
                             counts.only = FALSE,
                             verbose = TRUE, ...) {
    
    checkmate::assertClass(params, "splatPopParams")
    if (requireNamespace("vcfR", quietly = TRUE))
        
    if (verbose) {message("Getting parameters...")}
    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)
    
    
    sim.means <- splatPopSimulateMeans(vcf = vcf, 
                                       params = params,
                                       gff = gff,
                                       key = key,
                                       verbose = verbose)
    
    sim.sc <- splatPopSimulateSC(sim.means = sim.means$means, 
                                 params = params, 
                                 method = method,
                                 counts.only = counts.only, 
                                 verbose = verbose)
    
    out <- list(key = sim.means$key, means = sim.means$means, sc = sim.sc)
    
    return(out)
}


#' splatPopSimulateMeans
#'
#' Simulate mean expression levels for all genes for all samples, with between 
#' sample correlation structure simulated with eQTL effects and with the option
#' to simulate multiple groups (i.e. cell-types).
#'
#' @param vcf vcfR object containing genotypes of samples to simulate.
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details. 
#' @param gff Either NULL or a data.frame object containing a GFF/GTF file.
#' @param key Either FALSE or a data.frame object containing a full or partial
#'        splatPop key.
#' @param verbose logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details splatPopParams can be set in a variety of ways. 1. If
#' not provided, default parameters are used. 2. Default parameters can be
#' overridden by supplying desired parameters using \code{\link{setParams}}.
#' 3. Parameters can be estimated from real data of your choice using
#' \code{\link{splatPopEstimate}}.
#'
#' `splatPopSimulateMeans` involves the following steps:
#' \enumerate{
#'     \item Load population key or generate random or GFF/GTF based key.
#'     \item Format and subset genotype data from the VCF file.
#'     \item If not in key, assign expression mean and variance to each gene.
#'     \item If not in key, assign eGenes-eSNPs pairs and effect sizes.
#'     \item If not in key and groups >1, assign subset of eQTL associations as
#'     group-specific and assign DEG group effects.
#'     \item Simulate mean gene expression matrix without eQTL effects
#'     \item Quantile normalize by sample to fit single-cell expression
#'     distribution as defined in `splatEstimate`.
#'     \item Add quantile normalized gene mean and cv info the eQTL key.
#'     \item Add eQTL effects to means matrix.}
#'
#' @return A list containing: `means` a data.frame (or list of data.frames if
#' n.groups > 1) with the simulated mean gene expression value for each gene
#' (row) and each sample (column) and `key` a data.frame with population 
#' information including eQTL and group effects.
#'
#' @seealso
#' \code{\link{splatPopParseVCF}}, \code{\link{splatPopParseGenes}},
#' \code{\link{splatPopAssignMeans}},
#' \code{\link{splatPopQuantNorm}}, \code{\link{splatPopQuantNormKey}}
#' \code{\link{splatPopeQTLEffects}}, \code{\link{splatPopGroupEffects}}, 
#' \code{\link{splatPopSimMeans}}, \code{\link{splatPopSimEffects}}, 
#'
#' @export
#'

splatPopSimulateMeans <- function(vcf = mockVCF(), 
                                  params = newSplatPopParams(nGenes = 1000),
                                  verbose = TRUE, key = NULL, gff = NULL, ...){
    
    if (!requireNamespace("vcfR", quietly = TRUE)) {
        stop("The splatPop means simulation requires the 'vcfR' package.")}
    
    checkmate::assertClass(params, "splatPopParams")
    set.seed(getParam(params, "seed"))
    
    nGroups <- getParam(params, "nGroups")
    
    vcf.parsed <- splatPopParseVCF(vcf, params)
    group.names <- paste0("Group", seq_len(nGroups))

    # Genes from key or gff or mock (in that order)
    if(is.null(key)){key <- splatPopParseGenes(params, gff)}
    
    if (!all(c("exp_mean", "exp_cv") %in% names(key))){
        key <- splatPopAssignMeans(params, key)}
    
    if (!all(c("eQTL", "eSNP", "EffectSize") %in% names(key))){
        key <- splatPopeQTLEffects(params, key, vcf.parsed)
        
        if(length(group.names) > 1){
            key <- splatPopGroupEffects(params, key, group.names)
        }
    }
    
    if (verbose) {message("Simulating gene means for population...")}

    MeansPop <- splatPopSimMeans(vcf.parsed, key)

    MeansPop <- splatPopQuantNorm(params, MeansPop)
    key <- splatPopQuantNormKey(key, MeansPop) 

    eMeansPop <- splatPopSimEffects('global', key, vcf.parsed, MeansPop)

    if(length(group.names) > 1){
        eMeansPopq.groups <- list()
        for(id in group.names){
            eMeansPop.g <- splatPopSimEffects(id, key, vcf.parsed, eMeansPop)
            eMeansPop.g[eMeansPop.g <= 0] <- 1e-5 
            eMeansPopq.groups[[id]] <- eMeansPop.g}

        return(list(means = eMeansPopq.groups, key = key))
        
    }else{
        eMeansPop[eMeansPop <= 0] <- 1e-5 
        return(list(means = eMeansPop, key = key))
    }
}


#' splatPopSimulateSC
#'
#' Simulate count data for a population from a fictional single-cell 
#' RNA-seq experiment using the Splat method.
#' 
#' @param sim.means Matrix or list of matrices of gene means for the population. 
#'        Output from `splatPopSimulateMeans()`. 
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details. 
#' @param method which simulation method to use. Options are "single" which
#'        produces a single cell population for each sample, "groups" which 
#'        produces distinct groups (eg. cell types) for each sample (note, this
#'        creates separate groups from those created in `popSimulate` with only
#'        DE effects), and "paths" which selects cells from continuous
#'        trajectories (eg. differentiation processes).
#' @param counts.only logical. Whether to return only the counts.
#' @param verbose logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @return SingleCellExperiment with simulated data for whole population.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment cbind
#' @importFrom SummarizedExperiment rowData rowData<-
#' @export
#' 
splatPopSimulateSC <- function(sim.means,
                               params, 
                               method = c("single", "groups", "paths"),
                               counts.only = FALSE,
                               verbose = TRUE, ...){
    
    checkmate::assertClass(params, "splatPopParams")
    set.seed(getParam(params, "seed"))
    method <- match.arg(method)
    
    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)
    
    seed <- getParam(params, "seed")
    set.seed(seed)
    
    nGroups <- getParam(params, "nGroups")
    group.names <- paste0("Group", seq_len(nGroups))
    group.prop <- getParam(params, "group.prop")
    batchCells <- getParam(params, "batchCells")

    # Simulate sc counts with group-specific effects
    if (type(sim.means) == "list"){
        if(length(group.prop) != length(sim.means)){
            group.prop <- rep(1/length(sim.means), length(sim.means))}
        
        group.n <- lapply(group.prop, function(x) ceiling(x * batchCells))
        names(group.n) <- group.names
        samples <- names((sim.means[[1]]))
        group.sims <- list()
        
        for(g in group.names){
            if(verbose){message(paste0("Simulating sc counts for ", g, "..."))}
            paramsG <- params
            paramsG <- setParams(paramsG, batchCells = unlist(group.n[g]))
            
            sims <- lapply(samples, 
                           function(x) splatPopSimulateSample(params = paramsG, 
                                               method = method,
                                               sample.means = sim.means[[g]][x], 
                                               counts.only = counts.only,
                                               verbose = verbose))
            
            for(i in seq(1, length(sims))){
                s <- samples[i]
                sims[i][[1]]$Sample <- s
                sims[i][[1]]$Group <- g
                names(rowData(sims[i][[1]])) <- paste(s, g, 
                                                names(rowData(sims[i][[1]])), 
                                                sep = "_")}
            
            group.sims[[g]] <- do.call(SingleCellExperiment::cbind, sims)
        }
        
        sim.all <- do.call(SingleCellExperiment::cbind, group.sims)
        
    }else{
        if (verbose) {message("Simulating population single cell counts...")}
        samples <- names(sim.means)
        sims <- lapply(samples,
                       function(x) splatPopSimulateSample(params = params, 
                                                 method = method,
                                                 sample.means = sim.means[x], 
                                                 counts.only = counts.only,
                                                 verbose = verbose))
        for(i in seq(1, length(sims))){
            s <- samples[i]
            sims[i][[1]]$Sample <- s
            names(rowData(sims[i][[1]])) <- paste(s, 
                                                  names(rowData(sims[i][[1]])),
                                                  sep = "_")}
        
        sim.all <- do.call(SingleCellExperiment::cbind, sims)
    }
    
    sim.all <- splatPopCleanSCE(sim.all)

    if (verbose) {message("Done...")}
    return (sim.all)
}

#' splatPopSimulateSample simulation
#'
#' Simulate count data for one sample from a fictional single-cell RNA-seq 
#' experiment using the Splat method.
#'
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details.
#' @param method which simulation method to use. Options are "single" which
#'        produces a single population, "groups" which produces distinct groups
#'        (eg. cell types), "paths" which selects cells from continuous
#'        trajectories (eg. differentiation processes).
#' @param sample.means Gene means to use if running splatSimulatePop().
#' @param counts.only logical. Whether to return only the counts.
#' @param verbose logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' This function closely mirrors \code{\link{splatSimulate}}. The main
#' difference is that it takes the means simulated by splatPopSimulateMeans
#' instead of randomly sampling a mean for each gene. For details about this
#' function see the documentation for \code{\link{splatSimulate}}.
#'
#' @return SingleCellExperiment object containing the simulated counts and
#' intermediate values for one sample.
#'
#' @seealso
#' \code{\link{splatSimLibSizes}}, \code{\link{splatPopSimGeneMeans}},
#' \code{\link{splatSimBatchEffects}}, \code{\link{splatSimBatchCellMeans}},
#' \code{\link{splatSimDE}}, \code{\link{splatSimCellMeans}},
#' \code{\link{splatSimBCVMeans}}, \code{\link{splatSimTrueCounts}},
#' \code{\link{splatSimDropout}}, \code{\link{splatPopSimulateSC}}
#'
#' @importFrom SummarizedExperiment rowData colData colData<- assays
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods validObject
#' @export

splatPopSimulateSample <- function(params = newSplatPopParams(),
                                   method = c("single", "groups", "paths"),
                                   counts.only = FALSE,
                                   verbose = TRUE,
                                   sample.means, ...) {
    
    method <- match.arg(method)
    set.seed(getParam(params, "seed"))

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- nrow(sample.means)
    params <- setParams(params, nGenes = nGenes)
    nBatches <- getParam(params, "nBatches")
    batch.cells <- getParam(params, "batchCells")
    nGroups <- getParam(params, "nGroups")
    group.prob <- getParam(params, "group.prob")
    
    # Run sanity checks
    if (nGroups == 1 && method == "groups") {
        warning("nGroups is 1, switching to single mode")
        method <- "single"
    }

        # Set up name vectors
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- row.names(sample.means)
    batch.names <- paste0("Batch", seq_len(nBatches))
    if (method == "groups") {
        group.names <- paste0("Group", seq_len(nGroups))
    } else if (method == "paths") {
        group.names <- paste0("Path", seq_len(nGroups))
    }
    
    # Create SingleCellExperiment to store simulation
    cells <-  data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    features <- data.frame(Gene = gene.names)
    rownames(features) <- gene.names
    sim <- SingleCellExperiment(rowData = features, colData = cells,
                                metadata = list(Params = params))
    
    # Make batches vector which is the index of param$batchCells repeated
    # params$batchCells[index] times
    batches <- lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
                      b = batch.cells)
    batches <- unlist(batches)
    colData(sim)$Batch <- batch.names[batches]
    
    if (method != "single") {
        groups <- sample(seq_len(nGroups), nCells, prob = group.prob,
                         replace = TRUE)
        colData(sim)$Group <- factor(group.names[groups], levels = group.names)
    }
    
    sim <- splatSimLibSizes(sim, params)
    sim <- splatPopSimGeneMeans(sim, params, 
                                base.means.gene = sample.means[[1]])
    
    if (nBatches > 1) {sim <- splatSimBatchEffects(sim, params)}
    sim <- splatSimBatchCellMeans(sim, params)
    
    if (method == "single") {
        sim <- splatSimSingleCellMeans(sim, params)
    } else if (method == "groups") {
        sim <- splatSimGroupDE(sim, params)
        sim <- splatSimGroupCellMeans(sim, params)
    } else {
        sim <- splatSimPathDE(sim, params)
        sim <- splatSimPathCellMeans(sim, params)
    }
    
    sim <- splatSimBCVMeans(sim, params)
    sim <- splatSimTrueCounts(sim, params)
    sim <- splatSimDropout(sim, params)
    
    if (counts.only) {assays(sim)[!grepl('counts', names(assays(sim)))] <- NULL}
    
    return(sim)
    
}


#' Format and subset genotype data from a VCF file.
#'
#' Extract numeric alleles from vcf object and filter out SNPs missing genotype 
#' data or outside the Minor Allele Frequency range in `splatPopParams`.
#'
#' @param vcf vcfR object containing genotypes of samples to simulate.
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details.
#'
#' @return Genotype data.frame
#' 
#' @importFrom stats complete.cases na.omit
#' @importFrom utils data
#'
#' @export 

splatPopParseVCF <- function(vcf, params){
    
    if (requireNamespace("vcfR", quietly = TRUE)) 
    
    eqtl.maf.min <- getParam(params, "eqtl.maf.min")
    eqtl.maf.max <- getParam(params, "eqtl.maf.max")
    
    vcf@fix[,"ID"] <- NA
    vcf.gt <- as.data.frame(vcfR::extract.gt(vcf, 
                                             element = "GT",
                                             as.numeric = TRUE))
    vcf.gt$MAF <- maf(vcf)[,"Frequency"]
    
    # Filter SNPs with NAs and outside MAF range
    vcf.gt <- vcf.gt[complete.cases(vcf.gt), ]
    vcf.gt <- vcf.gt[vcf.gt$MAF >= eqtl.maf.min, ]
    vcf.gt <- vcf.gt[vcf.gt$MAF <= eqtl.maf.max, ]
    
    return(vcf.gt)
}

#' Generate population key matrix from random or gff provided gene information
#'
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details.
#' @param gff Either NULL or a data.frame object containing a GFF/GTF file.
#'
#' @return The Partial splatPop key data.frame.
#'
splatPopParseGenes <- function(params, gff){
    
    nGenes <- getParam(params, "nGenes")
    
    if(is.null(gff)){
        gff <- mockGFF(nGenes)
    }else{
        gff <- as.data.frame(gff)
        if ((length(names(gff)) < 8 | nrow(gff[gff[,3] == "gene",]) < 1)) {
            stop("GFF file did not contain gene features or other issue with
            file format. See example data.")}
        nGenes <- nrow(gff)
    }
    
    genes <- gff[gff[,3] == "gene",]
    genes$geneID <- paste0("gene_", formatC(1:nGenes, 
                                            width = nchar(nrow(genes)),
                                            format = "d",
                                            flag = "0"))
    genes$chr <- genes[, 1]
    genes$gene_start <- genes[, 4]
    genes$gene_end <- genes[, 5]
    genes$gene_mid <- floor(abs((genes$gene_start - genes$gene_end)/2)) + 
        genes$gene_start
    key <- genes[,c("geneID", "chr", "gene_start", "gene_end", "gene_mid")]
    
    return(key)
}

#' Sample expression mean and variance for each gene
#'
#' A mean and coefficient of variation is assigned to each gene by sampling from
#' gamma distributions parameterized from real data in `splatPopEstimate`.
#' The cv gamma distributions are binned by gene mean because the distribution
#' of variance in real data is not independent from the mean.
#'
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details.
#' @param key Partial splatPop key data.frame.
#'
#' @return The key updated with assigned means and variances.
#'
splatPopAssignMeans <- function(params, key){
    
    mean.shape <- getParam(params, "pop.mean.shape")
    mean.rate <- getParam(params, "pop.mean.rate")
    cv.param <- getParam(params, "pop.cv.param")
    
    # Sample gene means
    key$exp_mean <- rgamma(nrow(key), shape = mean.shape, rate = mean.rate)
    key$exp_cv <- NULL
    
    # Sample coefficient of variation for each gene
    for (g in 1:nrow(key)){
        exp.mean <- key[g, "exp_mean"]
        bin <- cv.param[(cv.param$start < exp.mean) &
                            (cv.param$end >= exp.mean), ]
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
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details.
#' @param key Partial splatPop key data.frame.
#' @param snps Data.frame of genotype information output from `vcf.parsed`.
#'
#' @return The key updated with assigned eQTL effects.
#'
splatPopeQTLEffects <- function(params, key, snps){
    
    eqtl.n <- getParam(params, "eqtl.n")
    if (eqtl.n > nrow(key)){eqtl.n <- nrow(key)} # Can't be greater than nGenes
    if (eqtl.n <= 1){eqtl.n <- nrow(key) * eqtl.n} # If <= 1 it is a percent
    
    eqtl.dist <- getParam(params, "eqtl.dist")
    eqtlES.shape <- getParam(params, "eqtl.ES.shape")
    eqtlES.rate <- getParam(params, "eqtl.ES.rate")
    
    # Set up data.frame to save info about selected eSNP-eGENE pairs
    snps.list <- row.names(snps)
    key.tmp <- key
    key$eQTL <- NA
    key$eSNP <- NA
    key$eSNP_MAF <- NA
    key$EffectSize <- 0
    
    for(i in seq_len(eqtl.n)){
        again <- TRUE
        while (again == TRUE){
            if(length(snps.list) == 0) {
                stop("Not enough SNPs in MAF range.")
            }
            s <- sample(snps.list, 1)
            snps.list <- snps.list[!snps.list == s]
            s.chr <- as.numeric(strsplit(s, "[_]")[[1]][1])
            s.loc <- as.numeric(strsplit(s, "[_]")[[1]][2])
            
            matches <- subset(key.tmp, (chr == s.chr & 
                                            gene_mid > s.loc - eqtl.dist &
                                            gene_mid < s.loc + eqtl.dist))
            if(nrow(matches) > 0){
                match <- sample(matches$geneID, 1)
                again <- FALSE
            }
        }
        
        key.tmp <- key.tmp[!(key.tmp$geneID == match),]
        ES <- rgamma(1, shape = eqtlES.shape, rate = eqtlES.rate)

        key[key$geneID == match, ]$eSNP <- s
        key[key$geneID == match, ]$EffectSize <- ES
        key[key$geneID == match, ]$eSNP_MAF <- snps[s, "MAF"]
        key[key$geneID == match, ]$eQTL <- "global"
        
        # Randomly make some effects negative
        key$EffectSize <- key$EffectSize * sample(c(1, -1),
                                                  length(key$EffectSize),
                                                  replace = TRUE)
    }
    
    return(key)
}


#' Assign group-specific eQTL and DEGs.
#'
#' If groups > 1, n eSNP-eGene pairs (n = 'eqtl.group.specific') are randomly
#' assigned as group specific.
#'
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details.
#' @param key Partial splatPop key data.frame.
#' @param groups array of group names
#'
#' @return The key updated with group eQTL and non-eQTL effects.
#'
splatPopGroupEffects <- function(params, key, groups){
    
    # Assign group-specific eQTL
    eqtl.n <- getParam(params, "eqtl.n")
    if (eqtl.n > nrow(key)){eqtl.n <- nrow(key)} # Can't be greater than nGenes
    if (eqtl.n <= 1){eqtl.n <- nrow(key) * eqtl.n} # If <= 1 it is a percent
    
    g.specific.perc <- getParam(params, "eqtl.group.specific")
    n.groups <- length(groups)
    n.specific.each <- ceiling(eqtl.n * g.specific.perc / n.groups)
    
    for(g in groups){
        glob.genes <- subset(key, eQTL == "global")$geneID
        g.specific <- sample(glob.genes, size = n.specific.each)
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
#' @param vcf vcfR object containing genotypes of samples to simulate.
#' @param key Partial splatPop key data.frame.
#'
#' @return Data.frame of gene mean expression levels WITHOUT eQTL effects.
#'
#' @importFrom stats rnorm
#'
splatPopSimMeans <- function(vcf, key){
    
    vcf$MAF <- NULL

    means <- lapply(key$geneID,
                    function(g) rnorm(ncol(vcf),
                                      mean = key[key$geneID == g,]$exp_mean,
                                      sd = key[key$geneID == g,]$exp_mean *
                                          key[key$geneID == g,]$exp_cv))
    
    means.df <- data.frame(do.call(rbind, means), row.names = key$geneID)
    names(means.df) <- names(vcf)
    
    return(means.df)
}


#' Add eQTL effects to means matrix
#'
#' Add eQTL effects and non-eQTL group effects to simulated means matrix.
#' The eQTL effects are incorporated using the following equation:
#' \deqn{Ygs = (ESg x Mgs x Gs) + Mgs }
#' Where Ygs is the mean for gene g and sample s, ESg is the effect size
#' assigned to g, Mgs is the mean expression assigned to g for s, and Gs
#' is the genotype (number of minor alleles) for s. Non-eQTL group effects are
#' incorporated as:
#' \deqn{Ygs = Mgs x GEg}
#' Where GEg is the group effect (i.e. differential expression) assigned to g.
#' To simulate multiple gene mean matrices with different group effects, this
#' function can be run with `id` designating the group id.
#'
#' @param id The group ID (e.g. "global" or "g1")
#' @param key Partial splatPop key data.frame.
#' @param snps Data.frame of genotype information output from `vcf.parsed`.
#' @param MeansPop Mean gene expression data.frame.
#'
#' @return Data.frame of gene mean expression levels WITH eQTL effects.
#'
splatPopSimEffects <- function(id, key, snps, MeansPop){
    
    # Add group-specific eQTL effects
    genes.use <- subset(key, eQTL == id)$geneID
    samples <- names(MeansPop)
    snps$MAF <- NULL
    
    for(g in genes.use){
        without.eqtl <- as.numeric(MeansPop[g,])
        ES <- key[key$geneID == g, "EffectSize"]
        eSNPsample <- key[key$geneID == g, "eSNP"]
        genotype <- as.numeric(snps[eSNPsample, samples])
        MeansPop[g,] <- without.eqtl + (ES * genotype * MeansPop[g,]) 
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
#' @param params splatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{splatPopParams}} for details.
#' @param MeansMatrix Mean gene expression with eQTL effects data.frame.
#'
#' @return Data.frame of quantile normalized gene mean expression levels.
#'
#' @export

splatPopQuantNorm <- function(params, MeansMatrix){
    
    if (requireNamespace("preprocessCore", quietly = TRUE))
    
    # Generate sample target distribution from sc parameters
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    target <- rgamma(10000, shape = mean.shape, rate = mean.rate)
    
    mat.norm <- preprocessCore::normalize.quantiles.use.target(
        as.matrix(MeansMatrix), target)
    mat.norm[mat.norm < 0] <- 0
    df.norm <- as.data.frame(mat.norm, row.names = row.names(MeansMatrix))
    names(df.norm) <- names(MeansMatrix)
    
    return(df.norm)
}

#' Add quantile normalized gene mean and cv info the eQTL key.
#'
#' @param key Partial splatPop key data.frame.
#' @param MeansMatrix The output from `splatPopQuantNorm` as a matrix or list of
#'                    matrices.
#'
#' @return Final eQTL key.
#'
#' @export
splatPopQuantNormKey <- function(key, MeansMatrix){
    
    if (type(MeansMatrix) == "list"){
        qn.means <- list()
        qn.cvs <- list()
        
        for(group in names(MeansMatrix)){
            qn.mean.gr <- rowMeans(MeansMatrix[[group]])
            qn.cv.gr <- apply(MeansMatrix[[group]], 1, FUN = co.var)
            qn.means[[group]] <- qn.mean.gr
            qn.cvs[[group]] <- qn.cv.gr
        }
        qn.mean <- rowMeans(as.data.frame(qn.means))
        qn.cv <- rowMeans(as.data.frame(qn.cvs))
        
    }else{
        qn.mean <- rowMeans(MeansMatrix)
        qn.cv <- apply(MeansMatrix, 1, FUN = co.var)
    }
    
    qn.df <- data.frame(list(expQN_mean = qn.mean, expQN_cv = qn.cv))
    qn.df$geneID <- row.names(qn.df)
    key <- merge(key, qn.df, by = "geneID")
    
    return(key)
}


#' Simulate gene means for splatPop
#'
#' Simulate outlier expression factors for splatPop. Genes with an outlier 
#' factor not equal to 1 are replaced with the median mean expression
#' multiplied by the outlier factor.
#'
#' @param sim SingleCellExperiment to add gene means to.
#' @param params SplatParams object with simulation parameters.
#' @param base.means.gene List of gene means for sample from matrix 
#'        generated by `splatPopSimulateMeans` and with the sample specified
#'        in `splatPopSimulateSC`.
#'
#' @return SingleCellExperiment with simulated gene means.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom stats rgamma median
splatPopSimGeneMeans <- function(sim, params, base.means.gene) {
    
    # Note: This function is similar to splatSimGeneMeans, except it uses the 
    # simulated gene mean instead of sampling one randomly. If changes are made
    # to the outlier method for splat, they should also be made here.
    
    nGenes <- getParam(params, "nGenes")
    out.prob <- getParam(params, "out.prob")
    out.facLoc <- getParam(params, "out.facLoc")
    out.facScale <- getParam(params, "out.facScale")

    # Add expression outliers
    outlier.facs <- getLNormFactors(nGenes, out.prob, 0, out.facLoc,
                                    out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]
    
    rowData(sim)$BaseGeneMean <- base.means.gene
    rowData(sim)$OutlierFactor <- outlier.facs
    rowData(sim)$GeneMean <- means.gene
    
    return(sim)
}

#' Clean up the population-scale SCE to remove redundant information
#'
#' @param sim.all SingleCellExperiment object with counts for all samples
#'
#' @return SingleCellExperiment with simulated sc counts.
#'
splatPopCleanSCE <- function(sim.all){
    
    # Remove redundant sce info
    keep.id <-  gsub("_Gene", "", names(rowData(sim.all))[1])
    
    shared.out.factor <- rowData(sim.all)[[paste0(keep.id, "_OutlierFactor")]]
    rowData(sim.all)[grepl("_OutlierFactor", names(rowData(sim.all)))] <- NULL
    rowData(sim.all)$Shared_OutlierFactor <- shared.out.factor
    
    rowData(sim.all)[grepl("_Gene", names(rowData(sim.all)))] <- NULL
    metadata(sim.all)[2:length(names(metadata(sim.all)))] <- NULL
    
    return(sim.all)
}
