#' @rdname newParams
#' @importFrom methods new
#' @export
newBASiCSParams <- function(...) {

    params <- new("BASiCSParams")
    params <- setParams(params, ...)

    return(params)
}

#' @importFrom checkmate checkInt checkDataFrame checkNumeric
setValidity("BASiCSParams", function(object) {

    v <- getParams(object, slotNames(object))

    nCells <- v$nCells
    nGenes <- v$nGenes
    checks <- c(nGenes = checkInt(v$nGenes, lower = 1),
                nCells = checkInt(v$nCells, lower = 1),
                seed = checkInt(v$seed, lower = 0),
                gene.params = checkDataFrame(v$gene.params,
                                             types = "numeric",
                                             any.missing = FALSE,
                                             min.rows = 1, ncols = 2),
                nSpikes = checkNumber(v$nSpikes, lower = 2, finite = TRUE),
                spike.means = checkNumeric(v$spike.means, lower = 0,
                                           finite = TRUE),
                cell.params = checkDataFrame(v$cell.params,
                                             types = "numeric",
                                             any.missing = FALSE,
                                             min.rows = 1, ncols = 2),
                theta = checkNumber(v$theta, lower = 0, finite = TRUE)
                )

    if (!all(colnames(v$gene.params) == c("Mean", "Delta"))) {
        checks <- c(checks, gene.params = "Incorrect column names")
    }

    if (!all(colnames(v$cell.params) == c("Phi", "S"))) {
        checks <- c(checks, cell.params = "Incorrect column names")
    }

    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

setMethod("show", "BASiCSParams", function(object) {

    pp <- list("Spike-ins:"   = c("[Number]" = "nSpikes",
                                  "[Means]"  = "spike.means"),
               "Variability:" = c("[Theta]"  = "theta"))

    callNextMethod()

    gene.params <- getParam(object, "gene.params")
    cell.params <- getParam(object, "cell.params")
    cat("Genes:", "\n")
    cat("(Params)", "\n")
    cat("data.frame with", dim(gene.params)[1], "features\n")
    print(head(gene.params, n = 3))
    cat("  ...  ...\n\n")

    cat("Cells:", "\n")
    cat("(Params)", "\n")
    cat("data.frame with", dim(cell.params)[1], "features\n")
    print(head(cell.params, n = 3))
    cat("  ...  ...\n\n")

    showPP(object, pp)
})
