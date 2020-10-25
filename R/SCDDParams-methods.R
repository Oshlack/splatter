#' @rdname newParams
#' @importFrom methods new
#' @export
newSCDDParams <- function(...) {

    checkDependencies("scDD")

    # Initialise scDatEx to avoid NOTE
    scDatEx <- NA
    utils::data("scDatEx", package = "scDD", envir = environment())

    params <- new("SCDDParams", SCdat = scDatEx)

    params <- setParams(params, ...)

    return(params)
}

#' @importFrom checkmate checkInt checkClass checkNumeric
setValidity("SCDDParams", function(object) {

    v <- getParams(object, slotNames(object))

    checks <- c(nGenes = checkInt(v$nGenes, lower = 1),
                nCells = checkInt(v$nCells, lower = 1),
                seed = checkInt(v$seed, lower = 0),
                SCDat = checkClass(v$SCdat, "SingleCellExperiment"),
                nDE = checkInt(v$nDE, lower = 0),
                nDP = checkInt(v$nDP, lower = 0),
                nDM = checkInt(v$nDM, lower = 0),
                nDB = checkInt(v$nDB, lower = 0),
                nEE = checkInt(v$nEE, lower = 0),
                nEP = checkInt(v$nEP, lower = 0),
                sd.range = checkNumeric(v$sd.range, lower = 0, len = 2),
                modeFC = checkNumeric(v$modeFC, lower = 0, len = 3),
                varInflation = checkNumeric(v$varInflation, lower = 0, len = 2),
                condition = checkCharacter(v$condition, len = 1))

    if (v$nGenes != (v$nDE + v$nDP + v$nDM + v$nDB + v$nEE + v$nEP)) {
        checks <- c(checks, nGenes = paste("nGenes is not consistent with",
                                           "nDE, nDP, nDM, nDB, nEE, nEP"))
    }

    if (!(v$condition %in% colnames(SummarizedExperiment::colData(v$SCdat)))) {
        checks <- c(checks, condition = paste("condition must be a column of",
                                              "colData(SCDat)"))
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
setMethod("setParam", "SCDDParams", function(object, name, value) {
    checkmate::assertString(name)

    if (name == "nGenes") {
        stop(name, " cannot be set directly, set nDE, nDP, nDM, nDB, nEE or ",
             "nEP instead")
    }

    nNames <- c("nDE", "nDP", "nDM", "nDB", "nEE", "nEP")
    if (name %in% nNames) {
        checkmate::assertInt(value, lower = 0)
        total <- value
        for (nName in nNames) {
            if (nName != name) {
                total <- total + getParam(object, nName)
            }
        }
        object <- setParamUnchecked(object, "nGenes", total)
    }

    object <- callNextMethod()

    return(object)
})

setMethod("show", "SCDDParams", function(object) {

    pp <- list("Genes:"       = c("(nDE)"       = "nDE",
                                  "(nDP)"       = "nDP",
                                  "(nDM)"       = "nDM",
                                  "(nDP)"       = "nDP",
                                  "(nEE)"       = "nEE",
                                  "(nEP)"       = "nEP"),
               "Fold change:" = c("[SD Range]"  = "sd.range",
                                  "[Mode FC]"   = "modeFC"),
               "Variance:"    = c("[Inflation]" = "varInflation"),
               "Condition:"   = c("[Condition]" = "condition"))

    callNextMethod()

    SCdat <- getParam(object, "SCdat")
    cat("Data:", "\n")
    cat("(SCdat)", "\n")
    cat("SingleCellExperiment with", dim(SCdat)[1], "features and",
        dim(SCdat)[2], "cells", "\n\n")

    showPP(object, pp)
})
