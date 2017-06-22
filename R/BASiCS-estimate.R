#' Estimate BASiCS simulation parameters
#'
#' Estimate simulation parameters for the BASiCS simulation from a real dataset.
#'
#' @param counts either a counts matrix or an SCESet object containing count
#'        data to estimate parameters from.
#' @param is.spike logical vector indicating which genes are spike-ins.
#' @param spike.input number of input molecules for each spike-in.
#' @param batch vector giving the batch that each cell belongs to.
#' @param n total number of MCMC iterations. Must be \code{>= max(4, thin)} and
#' a multiple of \code{thin}.
#' @param thin thining period for the MCMC sampler. Must be \code{>= 2}.
#' @param burn burn-in period for the MCMC sampler. Must be in the range
#' \code{1 <= burn < n} and a multiple of \code{thin}.
#' @param params BASiCSParams object to store estimated values in.
#' @param verbose logical. Whether to print progress messages.
#' @param progress logical. Whether to print additional BASiCS progress
#' messages.
#' @param ... Optional parameters passed to \code{\link[BASiCS]{BASiCS_MCMC}}.
#'
#' @details
#' This function is just a wrapper around \code{\link[BASiCS]{BASiCS_MCMC}} that
#' takes the output and converts it to a BASiCSParams object. See
#' \code{\link[BASiCS]{BASiCS_MCMC}} for details.
#'
#' @return BASiCSParams object containing the estimated parameters.
#'
#' @examples
#' data("sc_example_counts")
#' is.spike <- c(rep(FALSE, 40), rep(TRUE, 10))
#' spike.input <- rnorm(10, 500, 200)
#' params <- BASiCSEstimate(sc_example_counts[1:50, 1:20],
#'                          is.spike, spike.input)
#' params
#' @export
BASiCSEstimate <- function(counts, is.spike, spike.input,
                           batch = rep(1, ncol(counts)),
                           n = 40000, thin = 10, burn = 20000,
                           params = newBASiCSParams(), verbose = TRUE,
                           progress = TRUE, ...) {
    UseMethod("BASiCSEstimate")
}

#' @rdname BASiCSEstimate
#' @export
BASiCSEstimate.SCESet <- function(counts, is.spike, spike.input,
                                  batch = rep(1, ncol(counts)),
                                  n = 40000, thin = 10, burn = 20000,
                                  params = newBASiCSParams(), verbose = TRUE,
                                  progress = TRUE, ...) {
    counts <- scater::counts(counts)
    BASiCSEstimate(counts, params)
}

#' @rdname BASiCSEstimate
#' @export
BASiCSEstimate.matrix <- function(counts, is.spike, spike.input,
                                  batch = rep(1, ncol(counts)),
                                  n = 40000, thin = 10, burn = 20000,
                                  params = newBASiCSParams(), verbose = TRUE,
                                  progress = TRUE, ...) {

    checkmate::assertClass(params, "BASiCSParams")
    checkmate::assertNumeric(counts, lower = 0, finite = TRUE,
                             any.missing = FALSE)
    checkmate::assertLogical(is.spike, any.missing = FALSE, len = nrow(counts))
    checkmate::assertNumeric(spike.input, lower = 0, finite = TRUE,
                             len = sum(is.spike))
    checkmate::assertIntegerish(batch, lower = 0, any.missing = FALSE,
                                len = ncol(counts))
    checkmate::assertInt(thin, lower = 2)
    checkmate::assertInt(n, lower = max(4, thin))
    if ((n %% thin) != 0) {
        stop("'n' must be a multiple of 'thin'")
    }
    checkmate::assertInt(burn, lower = 1, upper = n - 1)
    if ((burn %% thin) != 0) {
        stop("'burn' must be a multiple of 'thin'")
    }

    spike.info <- data.frame(rownames(counts)[is.spike], spike.input)
    BASiCS.data <- suppressMessages(
                       BASiCS::newBASiCS_Data(counts, is.spike, spike.info,
                                              batch)
    )

    if (verbose) {
        mcmc <- BASiCS::BASiCS_MCMC(Data = BASiCS.data, N = n, Thin = thin,
                                    Burn = burn, PrintProgress = progress, ...)
    } else {
        mcmc <- suppressMessages(
                    BASiCS::BASiCS_MCMC(Data = BASiCS.data, N = n, Thin = thin,
                                        Burn = burn, PrintProgress = progress,
                                        ...)
        )
    }

    mcmc.summ <- BASiCS::Summary(mcmc)

    means <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "mu")[, 1]
    deltas <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "delta")[, 1]
    phis <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "phi")[, 1]
    ss <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "s")[, 1]
    thetas <- BASiCS::displaySummaryBASiCS(mcmc.summ, Param = "theta")[, 1]

    params <- setParams(params,
                        nGenes = sum(!is.spike),
                        batchCells = as.vector(table(batch)),
                        gene.params = data.frame(Mean = means, Delta = deltas),
                        nSpikes = sum(is.spike),
                        spike.means = spike.input,
                        cell.params = data.frame(Phi = phis, S = ss),
                        theta = thetas)

    return(params)
}
