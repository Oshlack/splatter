#' @rdname newParams
#' @importFrom methods new
#' @export
newBASiCSParams <- function(...) {

    checkDependencies("BASiCS")

    params <- new("BASiCSParams")
    params <- setParams(params, ...)

    return(params)
}

#' @importFrom checkmate checkInt checkDataFrame checkNumeric
setValidity("BASiCSParams", function(object) {

    object <- expandParams(object)
    v <- getParams(object, slotNames(object))

    nCells <- v$nCells
    nGenes <- v$nGenes
    nBatches <- v$nBatches
    checks <- c(seed = checkInt(v$seed, lower = 0),
                nGenes = checkInt(v$nGenes, lower = 1),
                nCells = checkInt(v$nCells, lower = 1),
                nBatches = checkInt(v$nBatches, lower = 1),
                batchCells = checkIntegerish(v$batchCells, lower = 1,
                                             len = nBatches),
                gene.params = checkDataFrame(v$gene.params,
                                             types = "numeric",
                                             any.missing = FALSE,
                                             min.rows = 1, ncols = 2),
                nSpikes = checkNumber(v$nSpikes, lower = 0, finite = TRUE),
                spike.means = checkNumeric(v$spike.means, lower = 0,
                                           finite = TRUE),
                cell.params = checkDataFrame(v$cell.params,
                                             types = "numeric",
                                             any.missing = FALSE,
                                             min.rows = 1, ncols = 2),
                theta = checkNumeric(v$theta, lower = 0, len = nBatches,
                                     finite = TRUE)
                )

    if (!all(colnames(v$gene.params) == c("Mean", "Delta"))) {
        checks <- c(checks, gene.params = "Incorrect column names")
    }

    if (!all(colnames(v$cell.params) == c("Phi", "S"))) {
        checks <- c(checks, cell.params = "Incorrect column names")
    }

    # Check batchCells matches nCells, nBatches
    if (nCells != sum(v$batchCells) || nBatches != length(v$batchCells)) {
        checks <- c(checks,
                    "nCells, nBatches and batchCells are not consistent")
    }

    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

#' @rdname setParam
setMethod("setParam", "BASiCSParams", function(object, name, value) {
    checkmate::assertString(name)

    if (name == "nCells" || name == "nBatches") {
        stop(name, " cannot be set directly, set batchCells instead")
    }

    if (name == "batchCells") {
        object <- setParamUnchecked(object, "nCells", sum(value))
        object <- setParamUnchecked(object, "nBatches", length(value))
    }

    object <- callNextMethod()

    return(object)
})

setMethod("show", "BASiCSParams", function(object) {

    pp <- list("Genes:"       = c("(Params)"      = "gene.params"),
               "Cells:"       = c("(Params)"      = "cell.params"),
               "Batches:"     = c("(Batches)"     = "nBatches",
                                  "(Batch Cells)" = "batchCells"),
               "Spike-ins:"   = c("(Number)"      = "nSpikes",
                                  "(Means)"       = "spike.means"),
               "Variability:" = c("(Theta)"       = "theta"))

    callNextMethod()

    showPP(object, pp)
})

#' @rdname expandParams
setMethod("expandParams", "BASiCSParams", function(object) {

    n <- getParam(object, "nBatches")

    vectors <- c("theta")

    object <- callNextMethod(object, vectors, n)

    return(object)
})
