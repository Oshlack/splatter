#' Estimate Lun2 simulation parameters
#'
#' Estimate simulation parameters for the Lun2 simulation from a real dataset.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param plates integer vector giving the plate that each cell originated from.
#' @param params Lun2Params object to store estimated values in.
#' @param min.size minimum size of clusters when identifying group of cells in
#'        the data.
#' @param verbose logical. Whether to show progress messages.
#' @param BPPARAM A \code{\link{BiocParallelParam}} instance giving the parallel
#'        back-end to be used. Default is \code{\link{SerialParam}} which uses a
#'        single core.
#'
#' @details
#' See \code{\link{Lun2Params}} for more details on the parameters.
#'
#' @return LunParams object containing the estimated parameters.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' library(scater)
#' set.seed(1)
#' sce <- mockSCE()
#'
#' plates <- as.numeric(factor(colData(sce)$Mutation_Status))
#' params <- lun2Estimate(sce, plates, min.size = 20)
#' params
#' }
#' @importFrom BiocParallel bplapply SerialParam
#' @export
lun2Estimate <- function(counts, plates, params = newLun2Params(),
                         min.size = 200, verbose = TRUE,
                         BPPARAM = SerialParam()) {
    UseMethod("lun2Estimate")
}

#' @rdname lun2Estimate
#' @export
lun2Estimate.SingleCellExperiment <- function(counts, plates,
                                              params = newLun2Params(),
                                              min.size = 200, verbose = TRUE,
                                              BPPARAM = SerialParam()) {
    counts <- getCounts(counts)
    lun2Estimate(counts, plates, params, min.size = min.size, verbose = verbose,
                 BPPARAM = BPPARAM)
}

#' @rdname lun2Estimate
#' @importFrom stats model.matrix
#' @importFrom locfit locfit
#' @export
lun2Estimate.matrix <- function(counts, plates, params = newLun2Params(),
                                min.size = 200, verbose = TRUE,
                                BPPARAM = SerialParam()) {

    checkDependencies("lun2")

    progress <- FALSE
    if (requireNamespace("progress", quietly = TRUE)) {
        progress <- TRUE
    }
    progress <- progress && verbose

    checkmate::assertClass(params, "Lun2Params")
    checkmate::assertInt(min.size, lower = 1, upper = length(plates))
    checkmate::assertIntegerish(plates, len = ncol(counts))

    if (length(unique(plates)) < 2) {
        stop("Plates must contain at least 2 values.")
    }

    dge <- edgeR::DGEList(counts)
    dge <- edgeR::`[.DGEList`(dge, rowMeans(dge$counts) >= 1, )

    # Estimate how many groups there are in the data
    if (verbose) {message("Estimating number of groups...")}
    groups <- scran::quickCluster(dge$counts, min.size)
    # Calculate normalisation factors
    if (verbose) {message("Computing normalisation factors...")}
    # Get the sizes for normalisation based on the number of cells in each
    # cluster
    min.cluster.size <- min(table(groups))
    if (min.cluster.size < 20) {
        sizes <- 10
    } else if (min.cluster.size < 100) {
        sizes <- seq(20, min.size, 20)
    } else {
        sizes <- seq(20, 100, 20)
    }
    sum.facs <- scran::computeSumFactors(dge$counts, cluster = groups,
                                         sizes = sizes, positive = TRUE)
    if (any(sum.facs == 0)) {
        warning("Some sum factors are zero. See ?scran::computeSumFactors ",
                "for details.")
        sum.facs <- sum.facs + 1e-6
    }
    dge$samples$norm.factors <- sum.facs / dge$samples$lib.size
    # Mean centre normalisation factors
    dge$samples$norm.factors <- dge$samples$norm.factors /
        exp(mean(log(dge$samples$norm.factors)))

    # Estimating the NB dispersion (assuming sufficient residual d.f. to
    # estimate the dispersion without EB shrinkage).
    if (verbose) {message("Estimating dispersions...")}
    plateX <- model.matrix(~plates)
    dge <- edgeR::estimateDisp(dge, plateX, prior.df = 0, trend = "none")

    # Estimating the log-overall mean
    if (verbose) {message("Estimating gene means...")}
    centered.off <- edgeR::getOffset(dge)
    centered.off <- centered.off - mean(centered.off)
    logmeans <- edgeR::mglmOneGroup(dge$counts, offset = centered.off,
                                    dispersion = dge$tagwise.dispersion)

    # Estimating the plate effect variance
    if (verbose) {message("Estimating plate effects...")}
    if (progress) {
        pb.format <- "[:bar] :percent eta: :eta"
        pb <- progress::progress_bar$new(format = pb.format,
                                         total = nrow(dge), clear = FALSE)
        pb$tick(0)
    } else if (verbose) {
        message("This may take some time. Install 'progress' to see a ",
                "progress bar.")
    }
    # As well as errors glmer produces warnings. Stop these showing because we
    # expect them.
    suppressWarnings(
    collected <- bplapply(seq_len(nrow(dge)), function(i) {
        if (progress) {pb$tick()}
        tryCatch({
            out <- lme4::glmer(
                       Counts ~ 0 + (1 | Plate) + offset(log(sum.facs)),
                       data = data.frame(Counts = as.integer(counts[i, ]),
                                         Group = groups,
                                         Plate = plates),
                       family = lme4::negative.binomial(1 / dge$tagwise[i]))
            output <- unlist(lme4::VarCorr(out))
            return(output)
        }, error = function(err) {
            output <- NA_real_
            return(output)
        })
    }, BPPARAM = BPPARAM))
    sigma2 <- mean(unlist(collected), na.rm = TRUE)

    # Repeating the estimation of the dispersion with ZINB models.
    if (verbose) {message("Estimating zero-inflated parameters...")}
    zinb.prop <- rep(-Inf, nrow(dge))
    zinb.disp <- dge$tagwise.dispersion
    zinb.mean <- exp(logmeans)
    nonzeros <- which(rowSums(dge$counts == 0) > 0)
    if (progress) {
        pb <- progress::progress_bar$new(format = "[:bar] :percent eta: :eta",
                                         total = length(nonzeros),
                                         clear = FALSE)
        pb$tick(0)
    } else if (verbose) {
        message("This may take some time. Install 'progress' to see a ",
                "progress bar.")
    }
    zinb.ests <- bplapply(nonzeros, function(i) {
        if (progress) {pb$tick()}
        zinb.est <- c(mean = zinb.mean[i], prop = zinb.prop[i],
                      disp = zinb.disp[i])
        tryCatch({
            zfit <- pscl::zeroinfl(dge$count[i, ] ~ 0 + plates | 1,
                                   dist = "negbin", offset = log(sum.facs))
            zinb.est <- c(mean = mean(exp(zfit$coefficients$count)),
                          prop = unname(zfit$coefficients$zero),
                          disp = 1 / zfit$theta)
        }, error = function(err) {})
        return(zinb.est)
    }, BPPARAM = BPPARAM)

    zinb.ests <- do.call("rbind", zinb.ests)

    zinb.prop[nonzeros] <- zinb.ests[, "prop"]
    zinb.disp[nonzeros] <- zinb.ests[, "disp"]
    zinb.mean[nonzeros] <- zinb.ests[, "mean"]

    zinb.prop <- exp(zinb.prop) / (1 + exp(zinb.prop))

    params <- setParams(params, nGenes = length(logmeans),
                        cell.plates = plates, plate.var = sigma2,
                        gene.params = data.frame(Mean = exp(logmeans),
                                                 Disp = dge$tagwise.dispersion),
                        zi.params = data.frame(Mean = zinb.mean,
                                               Disp = zinb.disp,
                                               Prop = zinb.prop),
                        cell.libSizes = dge$samples$lib.size)

    return(params)
}
