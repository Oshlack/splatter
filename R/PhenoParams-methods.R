#' @rdname newParams
#' @importFrom methods new
#' @export
newPhenoParams <- function(...) {

    checkDependencies("pheno")

    params <- new("PhenoParams")
    params <- setParams(params, ...)

    return(params)
}

setValidity("PhenoParams", function(object) {

    v <- getParams(object, slotNames(object))

    checks <- c(nGenes = checkmate::checkInt(v$nGenes, lower = 1),
                nCells = checkmate::checkInt(v$nCells, lower = 1),
                n.de = checkmate::checkInt(v$n.de, lower = 0),
                n.pst = checkmate::checkInt(v$n.pst, lower = 0),
                n.pst.beta = checkmate::checkInt(v$n.pst.beta, lower = 0),
                n.de.pst.beta = checkmate::checkInt(v$n.de.pst.beta, lower = 0),
                seed = checkmate::checkInt(v$seed, lower = 0))

    if (v$nGenes != (v$n.de + v$n.pst + v$n.pst.beta + v$n.de.pst.beta)) {
        checks <- c(checks,
                    nGenes = paste("nGenes is not consistent with",
                                   "n.de, n.pst, n.pst.beta, n.de.pst.beta"))
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
setMethod("setParam", "PhenoParams", function(object, name, value) {
    checkmate::assertString(name)

    if (name == "nGenes") {
        stop(name, " cannot be set directly, set n.de, n.pst, n.pst.beta or ",
             "n.de.pst.beta instead")
    }

    nNames <- c("n.de", "n.pst", "n.pst.beta", "n.de.pst.beta")
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

setMethod("show", "PhenoParams", function(object) {

    pp <- list("Genes:"  = c("[DE]"              = "n.de",
                             "[PST]"             = "n.pst",
                             "[PST + Beta]"      = "n.pst.beta",
                             "[DE + PST + Beta]" = "n.de.pst.beta"))

    callNextMethod()
    showPP(object, pp)
})
