#' @rdname newParams
#' @importFrom methods new
#' @export
newLun2Params <- function(...) {

    params <- new("Lun2Params")
    params <- setParams(params, ...)

    return(params)
}

#' @importFrom checkmate checkInt checkNumber checkNumeric checkDataFrame
#' checkCharacter checkFactor
setValidity("Lun2Params", function(object) {

    v <- getParams(object, slotNames(object))

    nCells <- v$nCells
    nGenes <- v$nGenes
    nPlates <- v$nPlates
    checks <- c(nGenes = checkInt(v$nGenes, lower = 1),
                nCells = checkInt(v$nCells, lower = 1),
                seed = checkInt(v$seed, lower = 0),
                nPlates = checkInt(v$nPlates, lower = 1),
                plate.ingroup = checkCharacter(v$plate.ingroup, min.len = 1),
                plate.mod = checkNumber(v$plate.mod, lower = 0),
                plate.var = checkNumber(v$plate.var, lower = 0),
                gene.params = checkDataFrame(v$gene.params,
                                             types = "numeric",
                                             any.missing = FALSE,
                                             min.rows = 1, ncols = 2),
                zi.params = checkDataFrame(v$zi.params,
                                           types = "numeric",
                                           any.missing = FALSE,
                                           min.rows = 1, ncols = 3),
                cell.plates = checkFactor(v$cell.plates, len = nCells),
                cell.libSizes = checkNumeric(v$cell.libSizes, lower = 0,
                                             min.len = 1, any.missing = FALSE,
                                             finite = TRUE),
                cell.libMod = checkNumber(v$cell.libMod, lower = 0),
                de.nGene = checkInt(v$de.nGenes, lower = 0),
                de.fc = checkNumber(v$de.fc, lower = 0))

    if (!all(colnames(v$gene.params) == c("Mean", "Disp"))) {
        checks <- c(checks, gene.params = "Incorrect column names")
    }

    if (!all(colnames(v$zi.params) == c("Mean", "Disp", "Prop"))) {
        checks <- c(checks, gene.params = "Incorrect column names")
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
#' @importFrom methods slotNames
setMethod("setParam", "Lun2Params", function(object, name, value) {
    checkmate::assertString(name)

    if (name == "nCells" || name == "nPlates") {
        stop(name, " cannot be set directly, set cell.plates instead")
    }

    if (name == "cell.plates") {
        object <- setParamUnchecked(object, "nCells", length(value))
        object <- setParamUnchecked(object, "nPlates", length(unique(value)))
        value <- factor(value)
    }

    object <- callNextMethod()

    return(object)
})

setMethod("show", "Lun2Params", function(object) {

    pp <- list("Genes:"     = c("(Params)"        = "gene.params",
                                "(ZI Params)"     = "zi.params"),
               "Plates:"    = c("[Number]"        = "nPlates",
                                "[Modifier]"      = "plate.mod",
                                "(Variance)"      = "plate.var"),
               "Cells:"     = c("[Plates]"        = "cell.plates",
                                "(Library Sizes)" = "cell.libSizes",
                                "[Lib Size Mod]"  = "cell.libMod"),
               "Diff Expr:" = c("[Genes]"         = "de.nGenes",
                                "[Fold change]"   = "de.fc"))

    callNextMethod()

    showPP(object, pp)
})
