#' Estimate Lun2 simulation parameters
#'
#' Estimate simulation parameters for the Lun2 simulation from a real dataset.
#'
#' @param counts either a counts matrix or an SCESet object containing count
#'        data to estimate parameters from.
#' @param plates integer vector giving the plate that each cell originated from.
#' @param params Lun2Params object to store estimated values in.
#' @param min.size minimum size of clusters when identifying group of cells in
#'        the data.
#' @param verbose logical. Whether to show progress messages.
#'
#' @details
#' See \code{\link{Lun2Params}} for more details on the parameters.
#'
#' @return LunParams object containing the estimated parameters.
#'
#' @examples
#' \dontrun{
#' data("sc_example_counts")
#' data("sc_example_cell_info")
#' plates <- factor(sc_example_cell_info$Mutation_Status)
#' params <- lun2Estimate(sc_example_counts, plates, min.size = 20)
#' params
#' }
#' @export
lun2Estimate <- function(counts, plates, params = newLun2Params(),
                         min.size = 200, verbose = TRUE) {
    UseMethod("lun2Estimate")
}

#' @rdname lun2Estimate
#' @export
lun2Estimate.SCESet <- function(counts, plates, params = newLun2Params(),
                                min.size = 200, verbose = TRUE) {
    counts <- scater::counts(counts)
    lun2Estimate(counts, plates, params, min.size = min.size, verbose = verbose)
}

#' @rdname lun2Estimate
#' @importFrom stats model.matrix
#' @export
lun2Estimate.matrix <- function(counts, plates, params = newLun2Params(),
                                min.size = 200, verbose = TRUE) {

    # Check suggested packages
    if (!requireNamespace("scran", quietly = TRUE)) {
        stop("The Lun2 simulation requires the 'scran' package for estimation.")
    }

    if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("The Lun2 simulation requires the 'lme4' package for estimation.")
    }

    if (!requireNamespace("pscl", quietly = TRUE)) {
        stop("The Lun2 simulation requires the 'pscl' package for estimation.")
    }

    progress <- FALSE
    if (requireNamespace("progress", quietly = TRUE)) {
        progress <- TRUE
    }
    progress <- progress && verbose

    checkmate::assertClass(params, "Lun2Params")
    checkmate::assertInt(min.size, lower = 20, upper = length(plates))
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
                                         sizes = sizes)
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
    collected <- lapply(seq_len(nrow(dge)), function(i) {
        if (progress) {pb$tick()}
        tryCatch({
            out <- lme4::glmer(Counts ~ 0 + (1 | Plate) + offset(log(sum.facs)),
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
    }))
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
    for (i in nonzeros) {
        if (progress) {pb$tick()}
        tryCatch({
            zfit <- pscl::zeroinfl(dge$count[i, ] ~ 0 + plates | 1,
                                   dist = "negbin", offset = log(sum.facs))
            zinb.mean[i] <- mean(exp(zfit$coefficients$count))
            zinb.prop[i] <- zfit$coefficients$zero
            zinb.disp[i] <- 1 / zfit$theta
        }, error = function(err) {})
    }
    zinb.prop <- exp(zinb.prop) / (1 + exp(zinb.prop))

    params <- setParams(params, nGenes = length(logmeans),
                        cell.plates = plates, plate.var = sigma2,
                        gene.means = exp(logmeans),
                        gene.disps = dge$tagwise.dispersion,
                        gene.ziMeans = zinb.mean, gene.ziDisps = zinb.disp,
                        gene.ziProps = zinb.prop,
                        cell.libSizes = dge$samples$lib.size)

    return(params)
}