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
    #    old.nCells <- getParam(object, "nCells")
        object <- setParamUnchecked(object, "nCells", length(value))
        object <- setParamUnchecked(object, "nPlates", length(unique(value)))
    #    if (length(value) != old.nCells) {
    #        warning("nCells has been changed. cell.libSizes will be sampled ",
    #                "to length nCells")
    #        selected <- sample(seq_len(old.nCells), length(value),
    #                           replace = TRUE)
    #        old.libSizes <- getParam(object, "cell.libSizes")
    #        object <- setParamUnchecked(object, "cell.libSizes",
    #                                    old.libSizes[selected])
    #    }
        value <- factor(value)
    }

    #if (name == "nGenes") {
    #    old.nGenes <- getParam(object, "nGenes")
    #    if (value != old.nGenes) {
    #        warning("nGenes has been changed. Gene parameter vectors will be ",
    #                "sampled to length new nGenes.")
    #        selected <- sample(seq_len(old.nGenes), size = value,
    #                           replace = TRUE)
    #        for (parameter in grep("gene", slotNames(object), value = TRUE)) {
    #            old.value <- getParam(object, parameter)
    #            object <- setParamUnchecked(object, parameter,
    #                                        old.value[selected])
    #        }
    #    }
    #}

    object <- callNextMethod()

    return(object)
})

setMethod("show", "Lun2Params", function(object) {

    pp <- list("Plates:"    = c("[Number]"        = "nPlates",
                                "[Modifier]"      = "plate.mod",
                                "(Variance)"      = "plate.var"),
               "Cells:"     = c("[Plates]"        = "cell.plates",
                                "(Library Sizes)" = "cell.libSizes",
                                "[Lib Size Mod]"  = "cell.libMod"),
               "Diff Expr:" = c("[Genes]"         = "de.nGenes",
                                "[Fold change]"   = "de.fc"))

    callNextMethod()

    gene.params <- getParam(object, "gene.params")
    zi.params <- getParam(object, "zi.params")
    cat("Genes:", "\n")
    cat("(Params)", "\n")
    cat("data.frame with", dim(gene.params)[1], "features\n")
    print(head(gene.params, n = 3))
    cat("  ...  ...\n")
    cat("(ZI Params)", "\n")
    cat("data.frame with", dim(zi.params)[1], "features\n")
    print(head(zi.params, n = 3))
    cat("  ...  ...  ...\n\n")

    showPP(object, pp)
})
