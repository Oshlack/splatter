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
#' @param eqtl Either NULL or if simulating population parameters directly from
#'        empirical data, a data.frame with empirical/desired eQTL results.
#'        To see required format, run `mockEmpiricalSet()` and see eqtl output.
#' @param means Either NULL or if simulating population parameters directly from
#'        empirical data, a Matrix of real gene means across a population, where
#'        each row is a gene and each column is an individual in the population.
#'        To see required format, run `mockEmpiricalSet()` and see means output.
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
#' \donttest{if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
#'     requireNamespace("preprocessCore", quietly = TRUE)) {
#'     vcf <- mockVCF()
#'     gff <- mockGFF()
#'     sim <- splatPopSimulate(vcf = vcf, gff = gff, sparsify = FALSE)
#' }}
#'
#' @export
splatPopSimulate <- function(params = newSplatPopParams(nGenes = 50),
                             vcf = mockVCF(),
                             method = c("single", "groups", "paths"),
                             gff = NULL,
                             eqtl = NULL,
                             means = NULL,
                             key = NULL,
                             counts.only = FALSE,
                             sparsify = TRUE,
                             verbose = TRUE, ...) {

    if (verbose) {message("Designing population...")}
    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)


    sim.means <- splatPopSimulateMeans(vcf = vcf,
                                       params = params,
                                       gff = gff,
                                       key = key,
                                       eqtl = eqtl,
                                       means = means,
                                       verbose = verbose)

    sim.sc <- splatPopSimulateSC(sim.means = sim.means$means,
                                 params = params,
                                 key = sim.means$key,
                                 conditions = sim.means$conditions,
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
#' @param eqtl Either NULL or if simulating population parameters directly from
#'        empirical data, a data.frame with empirical/desired eQTL results.
#'        To see required format, run `mockEmpiricalSet()` and see eqtl output.
#' @param means Either NULL or if simulating population parameters directly from
#'        empirical data, a Matrix of real gene means across a population, where
#'        each row is a gene and each column is an individual in the population.
#'        To see required format, run `mockEmpiricalSet()` and see means output.
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
#' (row) and each sample (column), `key` a data.frame with population
#' information including eQTL and group effects, and `condition` a named array
#' containing conditional group assignments for each sample.
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
                                  verbose = TRUE, key = NULL, gff = NULL,
                                  eqtl = NULL, means = NULL, ...){

    seed <- getParam(params, "seed")
    withr::with_seed(seed, {

    nGroups <- getParam(params, "nGroups")
    quant.norm <- getParam(params, "pop.quant.norm")

    vcf <- splatPopParseVCF(vcf, params)
    group.names <- paste0("Group", seq_len(nGroups))
    samples <- colnames(VariantAnnotation::geno(vcf)$GT)
    conditions <- splatPopDesignConditions(params, samples)

    # Genes from key if provided, or from empirical data or simulated from gff
    if (is.null(key)) {
        if (is.null(eqtl) || is.null(means)){
            if (verbose) {message("Simulating data for genes in GFF...")}
            key <- splatPopParseGenes(params, gff)
        }else {
            if (verbose) {message("Using base gene means from data provided...")}
            key <- splatPopParseEmpirical(vcf = vcf, gff = gff, eqtl = eqtl,
                                          means = means, params = params)
            params <- setParams(params, nGenes = nrow(key))
        }
    } else {
        if (verbose) {message("Simulating data for genes in key...")}
    }

    if (!all(c("meanSampled", "cvSampled") %in% names(key))) {
        key <- splatPopAssignMeans(params, key)
    }

    if (!all(c("eQTL.group", "eSNP.ID", "eQTL.EffectSize") %in% names(key))) {
        key <- splatPopeQTLEffects(params, key, vcf)
    }

    if (length(group.names) > 1) {
        key <- splatPopGroupEffects(params, key, group.names)
    }

    if (!all(c("eQTL.condition", "ConditionDE.Condition1") %in% names(key))) {
        key <- splatPopConditionEffects(params, key, conditions)
    }

    if (verbose) {message("Simulating gene means for population...")}

    means.pop <- splatPopSimMeans(vcf, key, means)

    if (quant.norm && ncol(means.pop) > 4){
        means.pop <- splatPopQuantNorm(params, means.pop)
        key <- splatPopQuantNormKey(key, means.pop)
    }

    eMeansPop <- splatPopSimEffects("global", key, conditions, vcf, means.pop)

    if (length(group.names) > 1) {
        eMeansPopq.groups <- list()
        for (id in group.names) {
            eMeansPop.g <- splatPopSimEffects(id, key, conditions, vcf,
                                              eMeansPop)
            eMeansPop.g[eMeansPop.g <= 0] <- 1e-5
            eMeansPopq.groups[[id]] <- eMeansPop.g
        }
        eMeansPop <- eMeansPopq.groups
    }

    sim.means <- splatPopSimConditionalEffects(key, eMeansPop, conditions)
    })

    return(list(means = sim.means, key = key, conditions=conditions))
}


#' splatPopParseEmpirical
#'
#' Parse splatPop key information from empirical data provided.
#'
#' NOTE: This function will cause some of the parameters in the splatPopParams
#' object to be ignored, such as population level gene mean and variance and
#' eQTL parameters.
#'
#' @param vcf VariantAnnotation object containing genotypes of samples.
#' @param gff Either NULL or a data.frame object containing a GFF/GTF file.
#' @param eqtl Either NULL or if simulating population parameters directly from
#'        empirical data, a data.frame with empirical/desired eQTL results.
#'        To see required format, run `mockEmpiricalSet()` and see eqtl output.
#' @param means Either NULL or if simulating population parameters directly from
#'        empirical data, a Matrix of real gene means across a population, where
#'        each row is a gene and each column is an individual in the population.
#'        To see required format, run `mockEmpiricalSet()` and see means output.
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#'
#' @details This function will ignore a number of parameters defined in
#' splatPopParams, instead pulling key information directly from provided VCF,
#' GFF, gene means, and eQTL mapping result data provided.
#'
#' @return A partial splatPop `key`
#'
#' @export
splatPopParseEmpirical <- function(vcf = vcf, gff = gff, eqtl = eqtl,
                                   means = means, params = params){

    key <- data.frame(chromosome = gff[, 1],
                      geneStart = gff[, 4],
                      geneEnd = gff[, 5],
                      geneMiddle = floor(abs((gff[, 4] - gff[, 5])/2)) +
                          gff[, 4],
                      meanSampled.noOutliers = rowMeans(means),
                      OutlierFactor = 1,
                      meanSampled = rowMeans(means),
                      cvSampled = apply(means, 1, co.var),
                      eQTL.group = ifelse(is.na(eqtl$snpID), NA, "global"),
                      eQTL.condition = "global",
                      eSNP.ID = eqtl$snpID,
                      eSNP.chromosome = eqtl$snpCHR,
                      eSNP.loc = eqtl$snpLOC,
                      eSNP.MAF = eqtl$snpMAF,
                      eQTL.EffectSize = eqtl$slope)
    row.names(key) <- eqtl$geneID


    return(key)
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
#' @param conditions named array with conditional group assignment for each
#'        sample. Output from `splatPopSimulateMeans()`.
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
                               conditions = NULL,
                               sparsify = TRUE,
                               verbose = TRUE, ...){

    method <- match.arg(method)
    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)

    seed <- getParam(params, "seed")
    withr::with_seed(seed, {

    nGroups <- getParam(params, "nGroups")
    group.names <- paste0("Group", seq_len(nGroups))
    group.prob <- getParam(params, "group.prob")
    nConditions <- getParam(params, "nConditions")
    condition.prob <- getParam(params, "condition.prob")

    batchCells <- getParam(params, "batchCells")

    if (!is.list(sim.means)) {
        sim.means <- list(Group1 = sim.means)
    }
    if (length(group.prob) != length(sim.means)) {
        group.prob <- rep(1 / length(sim.means), length(sim.means))
    }
    samples <- colnames((sim.means[[1]]))

    if (is.null(conditions)){
      conditions <- splatPopDesignConditions(params, samples)
    }
    batches <- splatPopDesignBatches(params, samples, verbose)

    # Simulate single-cell counts for each group/cell-type
    group.n <- lapply(group.prob, function(x) {ceiling(x * batchCells)})
    names(group.n) <- group.names
    group.sims <- list()
    for (g in group.names) {
        if (verbose) {message(paste0("Simulating sc counts for ", g, "..."))}

        paramsG <- setParams(params, batchCells = unlist(group.n[g]))
        sims <- lapply(samples, function(x) {
            splatPopSimulateSample(params = paramsG,
                                   method = method,
                                   sample.means = sim.means[[g]][, x],
                                   batch = batches[[x]],
                                   counts.only = counts.only,
                                   verbose = verbose)
        })


        for (i in seq(1, length(sims))) {
            s <- samples[i]
            c <- conditions[s]
            sims[i][[1]]$Sample <- s
            sims[i][[1]]$Group <- g
            sims[i][[1]]$Condition <- c
            names(rowData(sims[i][[1]])) <- paste(s, g,
                                            names(rowData(sims[i][[1]])),
                                            sep = "_")}
        group.sims[[g]] <- do.call(SingleCellExperiment::cbind, sims)
    }

    sim.all <- do.call(SingleCellExperiment::cbind, group.sims)
    sim.all <- splatPopCleanSCE(sim.all)

    metadata(sim.all)$Simulated_Means <- sim.means
    rowData(sim.all) <- merge(rowData(sim.all), key, by = "row.names")
    rowData(sim.all)$Row.names <- NULL

    if (sparsify) {
        if (verbose) {message("Sparsifying assays...")}
        assays(sim.all) <- sparsifyMatrices(assays(sim.all), auto = TRUE,
                                            verbose = verbose)
    }

    colnames(sim.all) <- paste(sim.all$Sample, sim.all$Cell, sep=":")
    })

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
#' @param batch Batch number.
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
                                   batch = "batch1",
                                   counts.only = FALSE,
                                   verbose = TRUE,
                                   sample.means, ...) {

    method <- match.arg(method)

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

    batch.list <- unlist(strsplit(batch, ","))
    batch.sims <- list()

    for (b in batch.list){
        # Get the parameters we are going to use
        if (isTRUE(getParam(params, "nCells.sample"))) {
            nGroups <- getParam(params, "nGroups")
            nCells.shape <- getParam(params, "nCells.shape")
            nCells.rate <- getParam(params, "nCells.rate")
            nCells <- rgamma(1, shape = nCells.shape, rate = nCells.rate)
            nCells <- ceiling(nCells / nGroups)
        } else {
            nCells <- getParam(params, "batchCells")[as.numeric(
                gsub("[^0-9.-]", "", b))]
        }
        # Set up name vectors
        cell.names <- paste0("Cell", seq_len(nCells))
        gene.names <- names(sample.means)
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
        colData(sim)$Batch <- b

        if (method != "single") {
            groups <- sample(seq_len(nGroups), nCells, prob = group.prob,
                             replace = TRUE)
            colData(sim)$Group <- factor(group.names[groups],
                                         levels = group.names)
        }

        sim <- splatSimLibSizes(sim, params)
        rowData(sim)$GeneMean <- sample.means
        if (nBatches > 1) {
            sim <- splatPopSimBatchEffects(sim, params)
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

        batch.sims[[b]] <- sim
    }

    batch.sims <- do.call(SingleCellExperiment::cbind, batch.sims)

    if (counts.only) {
        assays(batch.sims)[!grepl("counts", names(assays(batch.sims)))] <- NULL
    }

    return(batch.sims)

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

    vcf <- vcf[!is.na(SummarizedExperiment::rowRanges(vcf)$MAF)]
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

    if (is.null(gff)) {
        gff <- mockGFF(nGenes)
    } else {
        gff <- as.data.frame(gff)
        if ((length(names(gff)) < 8 | nrow(gff[gff[,3] == "gene",]) < 1)) {
            stop("GFF file did not contain gene features or other issue with ",
                 "file format. See example data.")
        }
        nGenes <- nrow(gff)
    }

    key <- gff[gff[,3] %in% c("gene", "Gene"),]

    row.names(key) <- paste0("gene_", formatC(seq_len(nGenes),
                                              width = nchar(nrow(key)),
                                              format = "d", flag = "0"))
    if (ncol(gff) == 9) {
        if (all(grepl("ENSG", gff[, 9]))) {
            gene_names <- gsub(";.*", "", gff[, 9])
            gene_names <- gsub(".*:", "", gene_names)
            row.names(key) <- gene_names
        }
    }

    key[['chromosome']] <- key[, 1]
    key[['geneStart']] <- key[, 4]
    key[['geneEnd']] <- key[, 5]
    key[['geneMiddle']] <- floor(abs((key$geneStart - key$geneEnd)/2)) +
        key$geneStart
    key <- key[,c("chromosome", "geneStart", "geneEnd", "geneMiddle")]

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
    out.prob <- getParam(params, "out.prob")
    out.facLoc <- getParam(params, "out.facLoc")
    out.facScale <- getParam(params, "out.facScale")

    # Scale the CV rate parameters
    cv.param <- getParam(params, "pop.cv.param")
    similarity.scale <- getParam(params, "similarity.scale")
    cv.param$rate <- cv.param$rate * similarity.scale

    # Sample gene means
    base.means.gene <- rgamma(nrow(key), shape = mean.shape, rate = mean.rate)

    # Add expression outliers
    outlier.facs <- getLNormFactors(nrow(key), out.prob, 0, out.facLoc,
                                    out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]

    key[["meanSampled.noOutliers"]] <- base.means.gene
    key[["OutlierFactor"]] <- outlier.facs
    key[["meanSampled"]] <- means.gene

    # Sample coefficient of variation for each gene
    key[["cvSampled"]] <- NULL
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
    if (eqtl.n > nrow(key)) {
        eqtl.n <- nrow(key) # Can't be greater than nGenes
    }
    if (eqtl.n <= 1) {
        eqtl.n <- nrow(key) * eqtl.n # If <= 1 it is a proportion
    }

    key_order <- row.names(key)
    eqtl.dist <- getParam(params, "eqtl.dist")
    eqtlES.shape <- getParam(params, "eqtl.ES.shape")
    eqtlES.rate <- getParam(params, "eqtl.ES.rate")
    eqtl.coreg <- getParam(params, "eqtl.coreg")

    # Set up data.frame to save info about selected eSNP-eGENE pairs
    snps.list <- row.names(vcf)
    key2 <- key
    key[, c("eQTL.group", "eQTL.condition", "eSNP.ID",
            "eSNP.chromosome", "eSNP.loc", "eSNP.MAF")] <- NA
    key[, "eQTL.EffectSize"] <- 0

    for (i in seq_len(eqtl.n)) {
        try <- TRUE
        while(try){
            g <- sample(row.names(key2), 1)

            g.rgn <- GenomicRanges::GRanges(
                seqnames = S4Vectors::Rle(key[g, "chromosome"]),
                ranges = IRanges::IRanges(key[g, "geneMiddle"] - eqtl.dist,
                                          key[g, "geneMiddle"] + eqtl.dist))

            vcf.subset <- IRanges::subsetByOverlaps(
                SummarizedExperiment::rowRanges(vcf), g.rgn)

            if (length(vcf.subset) == 0) {
                key2 <- key2[row.names(key2) != g, ]
                message(paste("No SNP within eqtl.dist limit for:", g))
            } else {
                try <- FALSE
            }
        }

        s <- sample(names(vcf.subset), 1)

        key2 <- key2[row.names(key2) != g, ]
        ES <- rgamma(1, shape = eqtlES.shape, rate = eqtlES.rate)

        key[g, ]$eSNP.ID <- s
        key[g, ]$eSNP.chromosome <- as.character(GenomeInfoDb::seqnames(vcf[s]))
        key[g, ]$eSNP.loc <- BiocGenerics::start(vcf[s])
        key[g, ]$eQTL.EffectSize <- ES
        key[g, ]$eSNP.MAF <- SummarizedExperiment::rowRanges(vcf[s])$MAF
        key[g, ]$eQTL.group <- "global"
        key[g, ]$eQTL.condition <- "global"
    }

    if (eqtl.coreg > 0) {
        coregulated.n <- ceiling(eqtl.coreg*eqtl.n)
        keyCoreg <- key[!is.na(key$eSNP.ID), ]
        keyCoreg <- keyCoreg[sample(1:eqtl.n, coregulated.n), ]
        keyCoreg <- keyCoreg[order(keyCoreg$chromosome, keyCoreg$geneMiddle), ]

        esnpKeep <- row.names(keyCoreg[seq(1,length(keyCoreg[,1]),2), ])
        esnpReplace <- row.names(keyCoreg[seq(2,length(keyCoreg[,1]),2), ])
        esnpKeep <- esnpKeep[1:length(esnpReplace)]
        key[row.names(key) %in% esnpReplace, grep("eSNP.", names(key))] <-
            key[row.names(key) %in% esnpKeep, grep("eSNP.", names(key))]
    }

    # Randomly make some effects negative
    # Note that all eQTL with the same eSNP will have the same sign
    esnps <- unique(key$eSNP.ID)
    signs <- sample(c(1, -1), length(esnps), replace = TRUE)
    names(signs) <- esnps
    key$sign <- signs[key$eSNP.ID]
    key$eQTL.EffectSize <- key$eQTL.EffectSize * key$sign
    key$sign <- NULL

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
#' @return The key updated with group eQTL and DE effects.
splatPopGroupEffects <- function(params, key, groups){

    # Assign group-specific eQTL effects
    eqtl.n <- getParam(params, "eqtl.n")
    if (eqtl.n > nrow(key)) {
        eqtl.n <- nrow(key) # Can't be greater than nGenes
    }
    if (eqtl.n <= 1) {
        eqtl.n <- nrow(key) * eqtl.n # If <= 1, it is a proportion
    }

    g.specific.perc <- getParam(params, "eqtl.group.specific")
    n.groups <- length(groups)
    n.specific.each <- ceiling(eqtl.n * g.specific.perc / n.groups)

    for (g in groups) {
        glob.genes <- row.names(subset(key, key$eQTL.group == "global"))
        g.specific <- sample(glob.genes, size = n.specific.each)
        key[["eQTL.group"]][row.names(key) %in% g.specific] <- g
    }

    # Assign group-specific DE effects
    nGenes <- nrow(key)
    de.prob <- getParam(params, "de.prob")
    de.downProb <- getParam(params, "de.downProb")
    de.facLoc <- getParam(params, "de.facLoc")
    de.facScale <- getParam(params, "de.facScale")

    for (idx in seq_len(n.groups)) {
        de.facs <- getLNormFactors(nGenes, de.prob, de.downProb,
                                   de.facLoc, de.facScale)
        key[, paste0("GroupDE.", groups[idx])] <- de.facs
    }

    return(key)
}


#' Assign Condition-specific eQTL and DEGs.
#'
#' If nConditions > 1, n eSNP-eGene pairs (n = 'eqtl.condition.specific') are
#' randomly assigned as condition specific.
#'
#' @param params SplatPopParams object containing parameters for population
#'        scale simulations. See \code{\link{SplatPopParams}} for details.
#' @param key Partial splatPop key data.frame.
#' @param conditions array of condition names
#'
#' @return The key updated with conditional eQTL and DE effects.
splatPopConditionEffects <- function(params, key, conditions){

    condition.names <- unique(conditions)
    if (length(condition.names) == 1) {
        key$eQTL.condition <- "global"
        key$ConditionDE.Condition1 <- 1
    } else {
        # Assign group-specific eQTL effects
        eqtl.n <- getParam(params, "eqtl.n")
        if (eqtl.n > nrow(key)) {
            eqtl.n <- nrow(key) # Can't be > nGenes
        }
        if (eqtl.n <= 1) {
            eqtl.n <- nrow(key) * eqtl.n # If <= 1 it is a proprotion
        }

        c.specific.perc <- getParam(params, "eqtl.condition.specific")
        n.conditions <- length(condition.names)
        n.specific.each <- ceiling(eqtl.n * c.specific.perc / n.conditions)

        for (c in condition.names) {
            glob.genes <- row.names(subset(key, key$eQTL.condition == "global"))
            c.specific <- sample(glob.genes, size = n.specific.each)
            key[["eQTL.condition"]][row.names(key) %in% c.specific] <- c
        }

        # Assign group-specific DE effects
        nGenes <- nrow(key)
        cde.prob <- getParam(params, "cde.prob")
        cde.downProb <- getParam(params, "cde.downProb")
        cde.facLoc <- getParam(params, "cde.facLoc")
        cde.facScale <- getParam(params, "cde.facScale")

        for (idx in seq_len(n.conditions)) {
            cde.facs <- getLNormFactors(nGenes, cde.prob, cde.downProb,
                                        cde.facLoc, cde.facScale)
            key[, paste0("ConditionDE.", condition.names[idx])] <- cde.facs
        }
    }

    return(key)
}

#' Simulate mean gene expression matrix without eQTL effects
#'
#' Gene mean expression levels are assigned to each gene for each pair randomly
#' from a normal distribution parameterized using the mean and cv assigned to
#' each gene in the key. If gene means matrix is provided, those will be used
#' instead.
#'
#' @param vcf VariantAnnotation object containing genotypes of samples.
#' @param key Partial splatPop key data.frame.
#' @param means Null or matrix of gene means to use
#'
#' @return matrix of gene mean expression levels WITHOUT eQTL effects.
#'
#' @importFrom stats rnorm
splatPopSimMeans <- function(vcf, key, means){

    if (is.null(means)) {
        means <- matrix(rnorm(ncol(vcf) * nrow(key), mean = key$meanSampled,
                              sd = key$meanSampled * key$cvSampled),
                        nrow = nrow(key), ncol = ncol(vcf))

        rownames(means) <- row.names(key)
        colnames(means) <- colnames(VariantAnnotation::geno(vcf)$GT)
    } else {
        means <- as.matrix(means)
    }

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
#' @param conditions array of condition assignments for each sample
#' @param vcf VariantAnnotation object containing genotypes of samples.
#' @param means.pop Population mean gene expression matrix
#'
#' @return data.frame of gene mean expression levels WITH eQTL effects.
#'
splatPopSimEffects <- function(id, key, conditions, vcf, means.pop){

    # Add group-specific eQTL effects
    genes.use <- row.names(subset(key, key$eQTL.group == id))
    samples <- colnames(means.pop)

    for (g in genes.use) {
        without.eqtl <- as.numeric(means.pop[g, ])
        ES <- key[g, "eQTL.EffectSize"]
        eSNPsample <- key[g, "eSNP.ID"]

        genotype_code <- as.array(VariantAnnotation::geno(
            vcf[eSNPsample, samples])$GT)
        genotype <- lengths(regmatches(genotype_code,
                                       gregexpr("1", genotype_code)))

        condition <- key[g, "eQTL.condition"]
        if (condition == "global") {
            cond <- rep(1, length(samples))
        } else {
            cond <- ifelse(conditions[samples] == condition, 1, 0)
        }

        means.pop[g, ] <- without.eqtl + (ES * genotype * means.pop[g, ] * cond)
    }

    # Add group-specific non-eQTL effects
    if (id != "global") {
        means.pop <- means.pop * key[, paste0("GroupDE.", id)]
    }

    means.pop[means.pop < 0] <- 0

    return(means.pop)
}


#' Add conditional DE effects to means matrix
#'
#' @param key Partial splatPop key data.frame.
#' @param conditions array of condition assignments for each sample
#' @param means.pop matrix or list of matrices with gene means.
#'
#' @return data.frame of gene mean expression levels WITH conditional DE
#' effects.
#'
splatPopSimConditionalEffects <- function(key, means.pop, conditions){

    # Add group-specific eQTL effects
    condition.list <- unique(conditions)

    for (c in condition.list) {
        c.samples <- names(conditions[conditions == c])
        c.de <- key[[paste0("ConditionDE.", c)]]

        if (is.list(means.pop)) {
            for (i in names(means.pop)) {
                means.pop[[i]][, c.samples] <- mapply("*",
                                                      means.pop[[i]][, c.samples],
                                                      c.de)
                means.pop[[i]][means.pop[[i]] < 0] <- 1e-5
            }
        } else {
            means.pop[, c.samples] <- mapply("*", means.pop[, c.samples], c.de)

            means.pop[means.pop <= 0] <- 1e-5
        }
    }

    return(means.pop)
}


#' Add conditional DE effects to means matrix
#'
#' @param id The group ID (e.g. "global" or "g1")
#' @param key Partial splatPop key data.frame.
#' @param vcf VariantAnnotation object containing genotypes of samples.
#' @param means.pop Population mean gene expression matrix
#'
#' @return data.frame of gene mean expression levels WITH eQTL effects.
#'
splatPopConditionalEffects <- function(id, key, vcf, means.pop){

    # Add group-specific eQTL effects
    genes.use <- row.names(subset(key, key$eQTL.group == id))
    samples <- colnames(means.pop)

    for (g in genes.use) {
        without.eqtl <- as.numeric(means.pop[g, ])
        ES <- key[g, "eQTL.EffectSize"]
        eSNPsample <- key[g, "eSNP.ID"]
        genotype_code <- as.array(VariantAnnotation::geno(
            vcf[eSNPsample, samples])$GT
        )
        genotype <- lengths(regmatches(genotype_code,
                                       gregexpr("1", genotype_code)))

        means.pop[g, ] <- without.eqtl + (ES * genotype * means.pop[g, ])
    }

    # Add group-specific non-eQTL effects
    if (id != "global") {
        means.pop <- means.pop * key[, paste0("GroupDE.", id)]
    }

    means.pop[means.pop < 0] <- 0
    eMeansPop[eMeansPop <= 0] <- 1e-5

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

    if (is.list(means)) {
        qn.means <- list()
        qn.cvs <- list()

        for (group in names(means)) {
            qn.mean.gr <- rowMeans(means[[group]])
            qn.cv.gr <- apply(means[[group]], 1, FUN = co.var)
            qn.means[[group]] <- qn.mean.gr
            qn.cvs[[group]] <- qn.cv.gr
        }
        qn.mean <- rowMeans(as.data.frame(qn.means))
        qn.cv <- rowMeans(as.data.frame(qn.cvs))

    } else {
        qn.mean <- rowMeans(means)
        qn.cv <- apply(means, 1, FUN = co.var)
    }

    qn.df <- data.frame(meanQuantileNorm = qn.mean, cvQuantileNorm = qn.cv)
    qn.df <- qn.df[row.names(key), ]
    key <- cbind(key, qn.df)

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

#' Simulate batch effects
#'
#' Simulate batch effects. Batch effect factors for each batch are produced
#' using \code{\link{getLNormFactors}} and these are added along with updated
#' means for each batch.
#'
#' @param sim SingleCellExperiment to add batch effects to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated batch effects.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
splatPopSimBatchEffects <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nBatches <- getParam(params, "nBatches")
    batch.facLoc <- getParam(params, "batch.facLoc")
    batch.facScale <- getParam(params, "batch.facScale")
    batch.rmEffect <- getParam(params, "batch.rmEffect")

    if (length(batch.facLoc) == 1) {
        batch.facLoc <- rep(batch.facLoc, nBatches)
    }
    if (length(batch.facScale) == 1) {
        batch.facScale <- rep(batch.facScale, nBatches)
    }

    batch <- unique(colData(sim)$Batch)
    batch.num <- as.numeric(gsub("[^0-9.-]", "", batch))

    batch.facs <- getLNormFactors(nGenes, 1, 0.5, batch.facLoc[batch.num],
                                  batch.facScale[batch.num])

    if (batch.rmEffect) {
        batch.facs <- rep(1, length(batch.facs))
    }

    rowData(sim)[[paste0("BatchFac", batch)]] <- batch.facs

    return(sim)
}

#' Set up pooled experimental design
#'
#' @param params SplatParams object with simulation parameters.
#' @param samples List of samples from vcf.
#' @param verbose logical. Whether to print progress messages.
#'
#' @return Vector with batch assignments for each sample.
#'
splatPopDesignBatches <- function(params, samples, verbose){

    nBatches <- getParam(params, "nBatches")
    batch.size <- getParam(params, "batch.size")

    if (nBatches == 1) {
        batches = rep("Batch1", length(samples))
        names(batches) <- samples
    } else {

        if (length(samples) < batch.size) {
            if (verbose) {message("Changing batch.size to ", length(samples))}
            batch.size <- length(samples)
        }

        if (length(samples) > nBatches * batch.size) {
            warning(
                "Not enough batches requested to include all samples! ",
                "Please increase nBatches or batch.size."
            )
        }

        xInt <- floor((nBatches * batch.size) / length(samples))
        xRem <- (nBatches * batch.size) %% length(samples)
        if (xRem > 0) {
            counts <- sort(table(c(rep(samples, xInt),
                                   sample(samples)[1: xRem])),
                           decreasing = TRUE)
        } else {
            counts <- table(rep(samples, xInt))
        }

        try.design <- TRUE
        while(try.design) {
            all.batches <- rep(paste0("Batch", 1:nBatches), batch.size)
            batches <- list()
            try.again <- FALSE
            for (s in names(counts)) {
                batches.remaining <- unique(all.batches)
                if (length(batches.remaining) >= counts[s]) {
                    b <- sample(batches.remaining, counts[s], replace = FALSE)
                    all.batches <- all.batches[-match(b, all.batches)]
                    batches[[s]] <- b
                } else {
                    try.again <- TRUE
                }
            }
            try.design <- try.again
        }
    }

    return(batches)
}

#' Set up designed experiments conditions
#'
#' @param params SplatParams object with simulation parameters.
#' @param samples List of samples from vcf.
#'
#' @return Vector with condition assignments for each sample.
#'
splatPopDesignConditions <- function(params, samples){

    condition.prob <- getParam(params, "condition.prob")

    if (length(condition.prob) == 1) {
        conditions <- rep("Condition1", length(samples))
        names(conditions) <- samples
    } else {
        ll <- ceiling(condition.prob * length(samples))
        conditions <- unlist(lapply(seq_along(ll), function(x) rep(paste0(
            "Condition", x), ll[[x]])))
        conditions <- sample(conditions, length(samples))
        names(conditions) <- samples
    }

    return(conditions)
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

    rowData(sim.all)[grepl("_Gene", names(rowData(sim.all)))] <- NULL
    metadata(sim.all)[2:length(names(metadata(sim.all)))] <- NULL

    return(sim.all)
}
