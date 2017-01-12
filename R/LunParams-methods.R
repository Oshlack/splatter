#' @rdname newParams
#' @importFrom methods new
#' @export
newLunParams <- function(...) {

    params <- new("LunParams")
    params <- setParams(params, ...)

    return(params)
}

#' @importFrom checkmate checkInt checkIntegerish checkNumber checkNumeric
setValidity("LunParams", function(object) {

    object <- expandParams(object)
    v <- getParams(object, slotNames(object))

    nGroups <- v$nGroups
    checks <- c(nGenes = checkInt(v$nGenes, lower = 1),
                nCells = checkInt(v$nCells, lower = 1),
                nGroups = checkInt(v$nGroups, lower = 1),
                groupCells = checkIntegerish(v$groupCells, lower = 1,
                                             len = nGroups),
                mean.rate = checkNumber(v$mean.rate, lower = 0),
                mean.shape = checkNumber(v$mean.shape, lower = 0),
                count.disp = checkNumber(v$count.disp, lower = 0),
                de.nGenes = checkIntegerish(v$de.nGenes, lower = 0,
                                            len = nGroups),
                de.upProp = checkNumeric(v$de.upProp, lower = 0, upper = 1,
                                         len = nGroups),
                de.upFC = checkNumeric(v$de.upFC, lower = 0,
                                       len = nGroups),
                de.downFC = checkNumeric(v$de.downFC, lower = 0,
                                         len = nGroups),
                seed = checkInt(v$seed, lower = 0))

    # Check groupCells matches nCells, nGroups
    if (v$nCells != sum(v$groupCells) || nGroups != length(v$groupCells)) {
        checks <- c(checks,
                    "nCells, nGroups and groupCells are not consistent")
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
setMethod("setParam", "LunParams", function(object, name, value) {
    checkmate::assertString(name)

    if (name == "nCells" || name == "nGroups") {
      stop(name, " cannot be set directly, set groupCells instead")
    }

    if (name == "groupCells") {
      object <- setParamUnchecked(object, "nCells", sum(value))
      object <- setParamUnchecked(object, "nGroups", length(value))
    }

    object <- callNextMethod()

    return(object)
})

setMethod("show", "LunParams", function(object) {

    pp <- list("Groups:"    = c("[Groups]"      = "nGroups",
                                "[Group Cells]" = "groupCells"),
               "Diff expr:" = c("[Genes]"       = "de.nGenes",
                                "[Up Prop]"     = "de.upProp",
                                "[Up FC]"       = "de.upFC",
                                "[Down FC]"     = "de.downFC"))

    callNextMethod()
    showPP(object, pp)
})

#' @rdname expandParams
setMethod("expandParams", "LunParams", function(object) {

    n <- getParam(object, "nGroups")

    vectors <- c("de.nGenes", "de.upProp", "de.upFC", "de.downFC")

    object <- callNextMethod(object, vectors, n)

    return(object)
})
