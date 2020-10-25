#' splatPop simulation
#'
#' Simulate scRNA-seq count data using the splat model for a population of
#' individuals with correlation structure.
#'
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#' @param vcf VariantAnnotation object containing genotypes of samples.
#' @param method which simulation method to use. Options are "single" which
#'        produces a single population, "groups" which produces distinct groups
#'        (eg. cell types), "paths" which selects cells from continuous
#'        trajectories (eg. differentiation processes).
#' @param gff Either NULL or a data.frame object containing a GFF/GTF file.
#' @param key Either NULL or a data.frame object containing a full or partial
#'        splatPop key.
#' @param counts.only logical. Whether to save only counts in sce object.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
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
#' @return SingleCellExperiment object containing simulated counts,
#' intermediate values like the gene means simulated in `splatPopSimulateMeans`,
#' and information about the differential expression and eQTL effects assigned
#' to each gene.
#'
#' @examples
#'
#' \donttest{
#' if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
#'     requireNamespace("preprocessCore", quietly = TRUE)) {
#'     sim <- splatPopSimulate()
#' }
#' }
#'
#' @export
splatPopSimulate <- function(params = newSplatPopParams(nGenes = 1000),
                             vcf = mockVCF(),
                             method = c("single", "groups", "paths"),
                             gff = NULL,
                             key = NULL,
                             counts.only = FALSE,
                             sparsify = TRUE,
                             verbose = TRUE, ...) {

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
                                 key = sim.means$key,
                                 method = method,
                                 counts.only = counts.only,
                                 sparsify = sparsify,
                                 verbose = verbose)

    return(sim.sc)
}


#' splatPopSimulateMeans
#'
#' Simulate mean expression levels for all genes for all samples, with between
#' sample correlation structure simulated with eQTL effects and with the option
#' to simulate multiple groups (i.e. cell-types).
#'
#' @param vcf VariantAnnotation object containing genotypes of samples.
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#' @param gff Either NULL or a data.frame object containing a GFF/GTF file.
#' @param key Either FALSE or a data.frame object containing a full or partial
#'        splatPop key.
#' @param verbose logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details SplatPopParams can be set in a variety of ways. 1. If
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
#'     \item Add eQTL effects to means matrix.
#' }
#'
#' @return A list containing: `means` a matrix (or list of matrices if
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
#' @examples
#'
#' \donttest{
#' if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
#'     requireNamespace("preprocessCore", quietly = TRUE)) {
#'     means <- splatPopSimulateMeans()
#' }
#' }
#'
#' @export
splatPopSimulateMeans <- function(vcf = mockVCF(),
                                  params = newSplatPopParams(nGenes = 1000),
                                  verbose = TRUE, key = NULL, gff = NULL, ...){

    set.seed(getParam(params, "seed"))

    nGroups <- getParam(params, "nGroups")

    vcf <- splatPopParseVCF(vcf, params)
    group.names <- paste0("Group", seq_len(nGroups))

    # Genes from key or gff or mock (in that order)
    if (is.null(key)) {
        key <- splatPopParseGenes(params, gff)
    }

    if (!all(c("meanSampled", "cvSampled") %in% names(key))) {
        key <- splatPopAssignMeans(params, key)
    }

    if (!all(c("eQTL.type", "eSNP.ID", "eQTL.EffectSize") %in% names(key))) {
        key <- splatPopeQTLEffects(params, key, vcf)

        if(length(group.names) > 1){
            key <- splatPopGroupEffects(params, key, group.names)
        }
    }

    if (verbose) {message("Simulating gene means for population...")}

    means.pop <- splatPopSimMeans(vcf, key)

    means.pop <- splatPopQuantNorm(params, means.pop)
    key <- splatPopQuantNormKey(key, means.pop)

    eMeansPop <- splatPopSimEffects("global", key, vcf, means.pop)

    if (length(group.names) > 1) {
        eMeansPopq.groups <- list()

        for (id in group.names) {
            eMeansPop.g <- splatPopSimEffects(id, key, vcf, eMeansPop)
            eMeansPop.g[eMeansPop.g <= 0] <- 1e-5
            eMeansPopq.groups[[id]] <- eMeansPop.g
        }

        return(list(means = eMeansPopq.groups, key = key))

    } else {
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
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#' @param key data.frame object containing a full or partial splatPop key.
#'        Output from `splatPopSimulateMeans()`.
#' @param method which simulation method to use. Options are "single" which
#'        produces a single cell population for each sample, "groups" which
#'        produces distinct groups (eg. cell types) for each sample (note, this
#'        creates separate groups from those created in `popSimulate` with only
#'        DE effects), and "paths" which selects cells from continuous
#'        trajectories (eg. differentiation processes).
#' @param counts.only logical. Whether to return only the counts.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @return SingleCellExperiment object containing simulated counts,
#' intermediate values like the gene means simulated in `splatPopSimulateMeans`,
#' and information about the differential expression and eQTL effects assigned
#' to each gene.
#'
#' @examples
#'
#' \donttest{
#' if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
#'     requireNamespace("preprocessCore", quietly = TRUE)) {
#'     params <- newSplatPopParams()
#'     sim.means <- splatPopSimulateMeans()
#'     sim <- splatPopSimulateSC(sim.means$means, params, sim.means$key)
#' }
#' }
#'
#' @importFrom SingleCellExperiment SingleCellExperiment cbind
#' @importFrom SummarizedExperiment rowData rowData<-
#' @export
splatPopSimulateSC <- function(sim.means,
                               params,
                               key,
                               method = c("single", "groups", "paths"),
                               counts.only = FALSE,
                               sparsify = TRUE,
                               verbose = TRUE, ...){

    set.seed(getParam(params, "seed"))
    method <- match.arg(method)

    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)

    seed <- getParam(params, "seed")
    set.seed(seed)

    nGroups <- getParam(params, "nGroups")
    group.names <- paste0("Group", seq_len(nGroups))
    group.prob <- getParam(params, "group.prob")
    batchCells <- getParam(params, "batchCells")

    # Simulate single-cell counts with group-specific effects
    if (is.list(sim.means)){
        if (length(group.prob) != length(sim.means)) {
            group.prob <- rep(1 / length(sim.means), length(sim.means))
        }

        group.n <- lapply(group.prob, function(x) {ceiling(x * batchCells)})
        names(group.n) <- group.names
        samples <- colnames((sim.means[[1]]))
        group.sims <- list()

        for (g in group.names) {
            if (verbose) {message(paste0("Simulating sc counts for ", g, "..."))}
            paramsG <- setParams(params, batchCells = unlist(group.n[g]))

            sims <- lapply(samples, function(x) {
                splatPopSimulateSample(params = paramsG,
                                       method = method,
                                       sample.means = sim.means[[g]][, x],
                                       counts.only = counts.only,
                                       verbose = verbose)
            })

            for (i in seq(1, length(sims))) {
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
        samples <- colnames(sim.means)
        sims <- lapply(samples, function(x) {
            splatPopSimulateSample(params = params,
                                   method = method,
                                   sample.means = sim.means[, x],
                                   counts.only = counts.only,
                                   verbose = verbose)
        })

        for (i in seq(1, length(sims))) {
            s <- samples[i]
            sims[i][[1]]$Sample <- s
            names(rowData(sims[i][[1]])) <- paste(s,
                                                  names(rowData(sims[i][[1]])),
                                                  sep = "_")
        }

        sim.all <- do.call(SingleCellExperiment::cbind, sims)
    }

    sim.all <- splatPopCleanSCE(sim.all)

    metadata(sim.all)$Simulated_Means <- sim.means
    rowData(sim.all) <- merge(rowData(sim.all), key,
                              by.x = "row.names", by.y = "geneID")

    if (sparsify) {
        if (verbose) {message("Sparsifying assays...")}
        assays(sim.all) <- sparsifyMatrices(assays(sim.all), auto = TRUE,
                                            verbose = verbose)
    }

    if (verbose) {message("Done!")}
    return (sim.all)
}

#' splatPopSimulateSample simulation
#'
#' Simulate count data for one sample from a fictional single-cell RNA-seq
#' experiment using the Splat method.
#'
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
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
splatPopSimulateSample <- function(params = newSplatPopParams(),
                                   method = c("single", "groups", "paths"),
                                   counts.only = FALSE,
                                   verbose = TRUE,
                                   sample.means, ...) {

    method <- match.arg(method)
    set.seed(getParam(params, "seed"))

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- length(sample.means)
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
    gene.names <- names(sample.means)
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
    sim <- splatPopSimGeneMeans(sim, params, base.means.gene = sample.means)

    if (nBatches > 1) {
        sim <- splatSimBatchEffects(sim, params)
    }
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

    if (counts.only) {
        assays(sim)[!grepl('counts', names(assays(sim)))] <- NULL
    }

    return(sim)

}


#' Format and subset genotype data from a VCF file.
#'
#' Extract numeric alleles from vcf object and filter out SNPs missing genotype
#' data or outside the Minor Allele Frequency range in `SplatPopParams`.
#'
#' @param vcf VariantAnnotation object containing genotypes of samples.
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#'
#' @return Genotype data.frame
#'
#' @importFrom stats complete.cases na.omit
#' @importFrom utils data
splatPopParseVCF <- function(vcf, params){

    # Filter SNPs with NAs and outside MAF range
    eqtl.maf.min <- getParam(params, "eqtl.maf.min")
    eqtl.maf.max <- getParam(params, "eqtl.maf.max")

    SummarizedExperiment::rowRanges(vcf)$MAF <-
        VariantAnnotation::snpSummary(vcf)$a1Freq

    vcf <- vcf[SummarizedExperiment::rowRanges(vcf)$MAF >= eqtl.maf.min &
                   SummarizedExperiment::rowRanges(vcf)$MAF <= eqtl.maf.max]

    return(vcf)
}

#' Generate population key matrix from random or gff provided gene information
#'
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#' @param gff Either NULL or a data.frame object containing a GFF/GTF file.
#'
#' @return The Partial splatPop key data.frame.
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

    key <- gff[gff[,3] %in% c("gene", "Gene"),]
    key$geneID <- paste0("gene_", formatC(seq_len(nGenes),
                                          width = nchar(nrow(key)),
                                          format = "d", flag = "0"))

    key[['chromosome']] <- key[, 1]
    key[['geneStart']] <- key[, 4]
    key[['geneEnd']] <- key[, 5]
    key[['geneMiddle']] <- floor(abs((key$geneStart - key$geneEnd)/2)) +
        key$geneStart
    key <- key[,c("geneID", "chromosome", "geneStart", "geneEnd", "geneMiddle")]

    return(key)
}

#' Sample expression mean and variance for each gene
#'
#' A mean and coefficient of variation is assigned to each gene by sampling from
#' gamma distributions parameterized from real data in `splatPopEstimate`.
#' The cv gamma distributions are binned by gene mean because the distribution
#' of variance in real data is not independent from the mean. The degree of
#' similarity between individuals can be further tuned using the
#' similarity.scale parameter in `SplatPopParams`.
#'
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#' @param key Partial splatPop key data.frame.
#'
#' @return The key updated with assigned means and variances.
splatPopAssignMeans <- function(params, key){

    mean.shape <- getParam(params, "pop.mean.shape")
    mean.rate <- getParam(params, "pop.mean.rate")
    cv.param <- getParam(params, "pop.cv.param")
    similarity.scale <- getParam(params, "similarity.scale")

    # Scale the CV rate parameters
    cv.param$rate <- cv.param$rate * similarity.scale

    # Sample gene means
    key[["meanSampled"]] <- rgamma(nrow(key), shape = mean.shape, rate = mean.rate)
    key[["cvSampled"]] <- NULL

    # Sample coefficient of variation for each gene
    for (g in seq_len(nrow(key))) {
        exp.mean <- key[g, "meanSampled"]
        bin <- cv.param[(cv.param$start < exp.mean) &
                            (cv.param$end >= exp.mean), ]
        key[g,"cvSampled"] <- rgamma(1, shape = bin$shape, rate = bin$rate)
    }

    return(key)
}


#' Assign eGenes-eSNPs pairs and effect sizes.
#'
#' Randomly pairs N genes (eGene) a SNP (eSNP) within the window size
#' (eqtl.dist) and assigns each pair an effect size sampled from a gamma
#' distribution parameterized using the effect sizes from a real eQTL study.
#'
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#' @param key Partial splatPop key data.frame.
#' @param vcf VariantAnnotation object containing genotypes of samples.
#'
#' @return The key updated with assigned eQTL effects.
splatPopeQTLEffects <- function(params, key, vcf){

    eqtl.n <- getParam(params, "eqtl.n")
    if (eqtl.n > nrow(key)){eqtl.n <- nrow(key)} # Can't be greater than nGenes
    if (eqtl.n <= 1){eqtl.n <- nrow(key) * eqtl.n} # If <= 1 it is a percent

    eqtl.dist <- getParam(params, "eqtl.dist")
    eqtlES.shape <- getParam(params, "eqtl.ES.shape")
    eqtlES.rate <- getParam(params, "eqtl.ES.rate")

    # Set up data.frame to save info about selected eSNP-eGENE pairs
    snps.list <- row.names(vcf)
    key2 <- key
    key[["eQTL.type"]] <- NA
    key[["eSNP.ID"]] <- NA
    key[["eSNP.chromosome"]] <- NA
    key[["eSNP.loc"]] <- NA
    key[["eSNP.MAF"]] <- NA
    key[["eQTL.EffectSize"]] <- 0

    for(i in seq_len(eqtl.n)){
        again <- TRUE
        while (again == TRUE){
            if(length(snps.list) == 0) {
                stop("Not enough SNPs in MAF range.")
            }
            s <- sample(snps.list, 1)
            snps.list <- snps.list[!snps.list == s]
            s.chr <- as.character(GenomeInfoDb::seqnames(vcf[s])@values)
            s.loc <- BiocGenerics::start(vcf[s])

            matches <- subset(key2, (as.character(key2$chromosome) == s.chr &
                                         key2$geneMiddle > s.loc - eqtl.dist &
                                         key2$geneMiddle < s.loc + eqtl.dist))
            if(nrow(matches) > 0){
                match <- sample(matches$geneID, 1)
                again <- FALSE
            }
        }

        key2 <- key2[!(key2$geneID == match),]
        ES <- rgamma(1, shape = eqtlES.shape, rate = eqtlES.rate)

        key[key$geneID == match, ]$eSNP.ID <- s
        key[key$geneID == match, ]$eSNP.chromosome <-
            as.character(GenomeInfoDb::seqnames(vcf[s]))
        key[key$geneID == match, ]$eSNP.loc <- BiocGenerics::start(vcf[s])
        key[key$geneID == match, ]$eQTL.EffectSize <- ES
        key[key$geneID == match, ]$eSNP.MAF <-
            SummarizedExperiment::rowRanges(vcf[s])$MAF
        key[key$geneID == match, ]$eQTL.type <- "global"

        # Randomly make some effects negative
        key$eQTL.EffectSize <- key$eQTL.EffectSize * sample(c(1, -1),
                                                  length(key$eQTL.EffectSize),
                                                  replace = TRUE)
    }

    return(key)
}


#' Assign group-specific eQTL and DEGs.
#'
#' If groups > 1, n eSNP-eGene pairs (n = 'eqtl.group.specific') are randomly
#' assigned as group specific.
#'
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#' @param key Partial splatPop key data.frame.
#' @param groups array of group names
#'
#' @return The key updated with group eQTL and non-eQTL effects.
splatPopGroupEffects <- function(params, key, groups){

    # Assign group-specific eQTL
    eqtl.n <- getParam(params, "eqtl.n")
    if (eqtl.n > nrow(key)){eqtl.n <- nrow(key)} # Can't be greater than nGenes
    if (eqtl.n <= 1){eqtl.n <- nrow(key) * eqtl.n} # If <= 1 it is a percent

    g.specific.perc <- getParam(params, "eqtl.group.specific")
    n.groups <- length(groups)
    n.specific.each <- ceiling(eqtl.n * g.specific.perc / n.groups)

    for(g in groups){
        glob.genes <- subset(key, key$eQTL.type == "global")$geneID
        g.specific <- sample(glob.genes, size = n.specific.each)
        key[["eQTL.type"]][key$geneID %in% g.specific] <- g
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
        key[, paste0(groups[idx], ".GroupEffect")] <- de.facs
    }
    return(key)
}

#' Simulate mean gene expression matrix without eQTL effects
#'
#' Gene mean expression levels are assigned to each gene for each pair randomly
#' from a normal distribution parameterized using the mean and cv assigned to
#' each gene in the key.
#'
#' @param vcf VariantAnnotation object containing genotypes of samples.
#' @param key Partial splatPop key data.frame.
#'
#' @return matrix of gene mean expression levels WITHOUT eQTL effects.
#'
#' @importFrom stats rnorm
splatPopSimMeans <- function(vcf, key){

    means <- matrix(rnorm(ncol(vcf) * nrow(key), mean = key$meanSampled,
                          sd = key$meanSampled * key$cvSampled),
                    nrow = nrow(key), ncol = ncol(vcf))

    rownames(means) <- key$geneID
    colnames(means) <- colnames(VariantAnnotation::geno(vcf)$GT)

    return(means)
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
#' @param vcf VariantAnnotation object containing genotypes of samples.
#' @param means.pop Population mean gene expression matrix
#'
#' @return data.frame of gene mean expression levels WITH eQTL effects.
#'
splatPopSimEffects <- function(id, key, vcf, means.pop){

    # Add group-specific eQTL effects
    genes.use <- subset(key, key$eQTL.type == id)$geneID
    samples <- colnames(means.pop)

    for (g in genes.use) {
        without.eqtl <- as.numeric(means.pop[g, ])
        ES <- key[key$geneID == g, "eQTL.EffectSize"]
        eSNPsample <- key[key$geneID == g, "eSNP.ID"]
        genotype_code <- as.array(VariantAnnotation::geno(
            vcf[eSNPsample, samples])$GT
        )
        genotype <- lengths(regmatches(genotype_code,
                                       gregexpr("1", genotype_code)))

        means.pop[g, ] <- without.eqtl + (ES * genotype * means.pop[g, ])
    }

    # Add group-specific non-eQTL effects
    if (id != "global") {
        means.pop <- means.pop * key[, paste0(id, ".GroupEffect")]
    }

    means.pop[means.pop < 0] <- 0

    return(means.pop)
}


#' Quantile normalize by sample to fit sc expression distribution.
#'
#' For each sample, expression values are quantile normalized (qgamma)
#' using the gamma distribution parameterized from splatEstimate(). This ensures
#' the simulated gene means reflect the distribution expected from a sc dataset
#' and not a bulk dataset.
#'
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#' @param means Mean gene expression matrix with eQTL effects.
#'
#' @return matrix of quantile normalized gene mean expression levels.
#'
#' @examples
#'
#' if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
#'     requireNamespace("preprocessCore", quietly = TRUE)) {
#'     bulk.means <- mockBulkMatrix(n.genes = 100, n.samples = 100)
#'     bulk.qnorm <- splatPopQuantNorm(newSplatPopParams(), bulk.means)
#' }
#'
#' @export
splatPopQuantNorm <- function(params, means){

    # Generate sample target distribution from sc parameters
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    target <- rgamma(10000, shape = mean.shape, rate = mean.rate)

    means.norm <- preprocessCore::normalize.quantiles.use.target(means, target)
    means.norm[means.norm < 0] <- 0

    rownames(means.norm) <- rownames(means)
    colnames(means.norm) <- colnames(means)

    return(means.norm)
}

#' Add quantile normalized gene mean and cv info the eQTL key.
#'
#' @param key Partial splatPop key data.frame.
#' @param means matrix or list of matrices containing means from
#'        `splatPopQuantNorm`
#'
#' @return Final eQTL key.
splatPopQuantNormKey <- function(key, means){

    if (is.list(means)){
        qn.means <- list()
        qn.cvs <- list()

        for(group in names(means)){
            qn.mean.gr <- rowMeans(means[[group]])
            qn.cv.gr <- apply(means[[group]], 1, FUN = co.var)
            qn.means[[group]] <- qn.mean.gr
            qn.cvs[[group]] <- qn.cv.gr
        }
        qn.mean <- rowMeans(as.data.frame(qn.means))
        qn.cv <- rowMeans(as.data.frame(qn.cvs))

    }else{
        qn.mean <- rowMeans(means)
        qn.cv <- apply(means, 1, FUN = co.var)
    }

    qn.df <- data.frame(list(meanQuantileNorm = qn.mean,
                             cvQuantileNorm = qn.cv))
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
