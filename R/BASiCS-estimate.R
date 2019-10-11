#' Estimate BASiCS simulation parameters
#'
#' Estimate simulation parameters for the BASiCS simulation from a real dataset.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param spike.info data.frame describing spike-ins with two columns: "Name"
#'        giving the names of the spike-in features (must match
#'        \code{rownames(counts)}) and "Input" giving the number of input
#'        molecules.
#' @param batch vector giving the batch that each cell belongs to.
#' @param n total number of MCMC iterations. Must be \code{>= max(4, thin)} and
#' a multiple of \code{thin}.
#' @param thin thining period for the MCMC sampler. Must be \code{>= 2}.
#' @param burn burn-in period for the MCMC sampler. Must be in the range
#' \code{1 <= burn < n} and a multiple of \code{thin}.
#' @param regression logical. Whether to use regression to identify
#' over-dispersion. See \code{\link[BASiCS]{BASiCS_MCMC}} for details.
#' @param params BASiCSParams object to store estimated values in.
#' @param verbose logical. Whether to print progress messages.
#' @param progress logical. Whether to print additional BASiCS progress
#' messages.
#' @param ... Optional parameters passed to \code{\link[BASiCS]{BASiCS_MCMC}}.
#'
#' @details
#' This function is just a wrapper around \code{\link[BASiCS]{BASiCS_MCMC}} that
#' takes the output and converts it to a BASiCSParams object. Either a set of
#' spike-ins or batch information (or both) must be supplied. If only batch
#' information is provided there must be at least two batches. See
#' \code{\link[BASiCS]{BASiCS_MCMC}} for details.
#'
#' @return BASiCSParams object containing the estimated parameters.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' library(scater)
#' set.seed(1)
#' sce <- mockSCE()
#'
#' spike.info <- data.frame(Name = rownames(sce)[1:10],
#'                          Input = rnorm(10, 500, 200),
#'                          stringsAsFactors = FALSE)
#' params <- BASiCSEstimate(sce[1:100, 1:30], spike.info)
#' params
#' }
#' @export
BASiCSEstimate <- function(counts, spike.info = NULL, batch = NULL,
                           n = 20000, thin = 10, burn = 5000,
                           regression = TRUE,
                           params = newBASiCSParams(), verbose = TRUE,
                           progress = TRUE, ...) {
    UseMethod("BASiCSEstimate")
}

#' @rdname BASiCSEstimate
#' @export
BASiCSEstimate.SingleCellExperiment <- function(counts, spike.info = NULL,
                                                batch = NULL, n = 20000,
                                                thin = 10, burn = 5000,
                                                regression = TRUE,
                                                params = newBASiCSParams(),
                                                verbose = TRUE, progress = TRUE,
                                                ...) {
    counts <- getCounts(counts)
    BASiCSEstimate(counts, spike.info, batch, n, thin, burn, regression,
                   params, verbose, progress, ...)
}

#' @rdname BASiCSEstimate
#' @export
BASiCSEstimate.matrix <- function(counts, spike.info = NULL, batch = NULL,
                                  n = 20000, thin = 10, burn = 5000,
                                  regression = TRUE,
                                  params = newBASiCSParams(), verbose = TRUE,
                                  progress = TRUE, ...) {

    checkmate::assertClass(params, "BASiCSParams")
    checkmate::assertMatrix(counts, mode = "numeric", any.missing = FALSE,
                            min.rows = 1, min.cols = 1, row.names = "unique",
                            col.names = "unique")
    if (is.null(spike.info) && is.null(batch)) {
        stop("At least one of spike.info and batch must be provided")
    }
    if (!is.null(spike.info)) {
        checkmate::assertDataFrame(spike.info, any.missing = FALSE,
                                   min.rows = 1, ncols = 2)
        if (!all(colnames(spike.info) == c("Name", "Input"))) {
            stop("spike.info must have columns named 'Name' and 'Input'")
        }
        checkmate::assertCharacter(spike.info$Name, min.chars = 1,
                                   unique = TRUE)
        checkmate::assertNumeric(spike.info$Input, lower = 0, finite = TRUE)
    } #else {
    #    spike.info <- data.frame(Name = c(), Input = c())
    #}
    if (!is.null(batch)) {
        checkmate::assertIntegerish(batch, lower = 0, any.missing = FALSE,
                                    len = ncol(counts))
        if (is.null(spike.info) && length(unique(batch)) == 1) {
            stop("If spike.info is not provided there must be at least two ",
                 "batches")
        }
    } else {
        batch <- rep(1, ncol(counts))
    }
    checkmate::assertInt(thin, lower = 2)
    checkmate::assertInt(n, lower = max(4, thin))
    if ((n %% thin) != 0) {
        stop("'n' must be a multiple of 'thin'")
    }
    checkmate::assertInt(burn, lower = 1, upper = n - 1)
    if ((burn %% thin) != 0) {
        stop("'burn' must be a multiple of 'thin'")
    }
    checkmate::assertFlag(regression)

    is.spike <- rownames(counts) %in% spike.info$Name

    BASiCS.data <- suppressMessages(
                       BASiCS::newBASiCS_Data(counts, is.spike, spike.info,
                                              batch)
    )

    with.spikes <- sum(is.spike) > 1

    if (verbose) {
        mcmc <- BASiCS::BASiCS_MCMC(Data = BASiCS.data, N = n, Thin = thin,
                                    Burn = burn, Regression = regression,
                                    PrintProgress = progress,
                                    WithSpikes = with.spikes, ...)
    } else {
        mcmc <- suppressMessages(
                    BASiCS::BASiCS_MCMC(Data = BASiCS.data, N = n, Thin = thin,
                                        Burn = burn, Regression = regression,
                                        PrintProgress = progress,
                                        WithSpikes = with.spikes, ...)
        )
    }

    mcmc.summ <- BASiCS::Summary(mcmc)

    means <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "mu")[, 1]
    deltas <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "delta")[, 1]
    if (!is.null(spike.info)) {
        phis <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "phi")[, 1]
    } else {
        phis <- rep(1, ncol(counts))
    }
    ss <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "s")[, 1]
    thetas <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "theta")[, 1]

    params <- setParams(params,
                        nGenes = sum(!is.spike),
                        batchCells = as.vector(table(batch)),
                        gene.params = data.frame(Mean = means, Delta = deltas),
                        nSpikes = sum(is.spike),
                        spike.means = ifelse(!is.null(spike.info),
                                             spike.info$Input, numeric()),
                        cell.params = data.frame(Phi = phis, S = ss),
                        theta = thetas)

    return(params)
}
