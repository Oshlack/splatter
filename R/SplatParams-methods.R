#' @rdname newParams
#' @importFrom methods new
#' @export
newSplatParams <- function(...) {

    params <- new("SplatParams")
    params <- setParams(params, ...)

    return(params)
}

#' @importFrom checkmate checkInt checkIntegerish checkNumber checkNumeric
#' checkFlag
setValidity("SplatParams", function(object) {

    object <- expandParams(object)
    v <- getParams(object, c(slotNames(object)))

    nGroups <- v$nGroups
    checks <- c(nGenes = checkInt(v$nGenes, lower = 1),
                nCells = checkInt(v$nCells, lower = 1),
                nGroups = checkInt(v$nGroups, lower = 1),
                groupCells = checkIntegerish(v$groupCells, lower = 1,
                                             len = nGroups),
                mean.rate = checkNumber(v$mean.rate, lower = 0),
                mean.shape = checkNumber(v$mean.shape, lower = 0),
                lib.loc = checkNumber(v$lib.loc),
                lib.scale = checkNumber(v$lib.scale, lower = 0),
                out.prob = checkNumber(v$out.prob, lower = 0, upper = 1),
                out.facLoc = checkNumber(v$out.facLoc),
                out.facScale = checkNumber(v$out.facScale, lower = 0),
                de.prob = checkNumeric(v$de.prob, lower = 0, upper = 1,
                                       len = nGroups),
                de.downProb = checkNumeric(v$de.downProb, lower = 0, upper = 1,
                                           len = nGroups),
                de.facLoc = checkNumeric(v$de.facLoc, len = nGroups),
                de.facScale = checkNumeric(v$de.facScale, lower = 0,
                                           len = nGroups),
                bcv.common = checkNumber(v$bcv.common, lower = 0),
                bcv.df = checkNumber(v$bcv.df, lower = 0),
                dropout.present = checkFlag(v$dropout.present),
                dropout.mid = checkNumber(v$dropout.mid),
                dropout.shape = checkNumber(v$dropout.shape),
                path.from = checkIntegerish(v$path.from, lower = 0,
                                            upper = nGroups, len = nGroups),
                path.length = checkIntegerish(v$path.length, lower = 1,
                                               len = nGroups),
                path.skew = checkNumeric(v$path.skew, lower = 0, upper = 1,
                                         len = nGroups),
                path.nonlinearProb = checkNumber(v$path.nonlinearProb,
                                                 lower = 0, upper = 1),
                path.sigmaFac = checkNumber(v$path.sigmaFac, lower = 0),
                seed = checkInt(v$seed, lower = 0))

    # Check groupCells matches nCells, nGroups
    if (v$nCells != sum(v$groupCells) || nGroups != length(v$groupCells)) {
        checks <- c(checks,
                    "nCells, nGroups and groupCells are not consistent")
    }

    # Check path.from
    if (!(0 %in% v$path.from)) {
       checks <- c(checks, path.from = "origin must be specified in path.from")
    } else if (any(v$path.from == seq_len(nGroups))) {
        checks <- c(checks, stop("path cannot begin at itself"))
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
setMethod("setParam", "SplatParams",function(object, name, value) {
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

#' @importFrom methods callNextMethod
setMethod("show", "SplatParams", function(object) {

    pp <- list("Groups:"         = c("[Groups]"       = "nGroups",
                                     "[Group Cells]"  = "groupCells"),
               "Mean:"           = c("(Rate)"         = "mean.rate",
                                     "(Shape)"        = "mean.shape"),
               "Library size:"   = c("(Location)"     = "lib.loc",
                                     "(Scale)"        = "lib.scale"),
               "Exprs outliers:" = c("(Probability)"  = "out.prob",
                                     "(Location)"     = "out.facLoc",
                                     "(Scale)"        = "out.facScale"),
               "Diff expr:"      = c("[Probability]"  = "de.prob",
                                     "[Down Prob]"    = "de.downProb",
                                     "[Location]"     = "de.facLoc",
                                     "[Scale]"        = "de.facScale"),
               "BCV:"            = c("(Common Disp)"  = "bcv.common",
                                     "(DoF)"          = "bcv.df"),
               "Dropout:"        = c("(Present)"      = "dropout.present",
                                     "(Midpoint)"     = "dropout.mid",
                                     "(Shape)"        = "dropout.shape"),
               "Paths:"          = c("[From]"         = "path.from",
                                     "[Length]"       = "path.length",
                                     "[Skew]"         = "path.skew",
                                     "[Non-linear]"   = "path.nonlinearProb",
                                     "[Sigma Factor]" = "path.sigmaFac"))

    callNextMethod()
    showPP(object, pp)
})

#' @rdname expandParams
setMethod("expandParams", "SplatParams", function(object) {

    n <- getParam(object, "nGroups")

    vectors <- c("de.prob", "de.downProb", "de.facLoc", "de.facScale",
                 "path.from", "path.length", "path.skew")

    object <- callNextMethod(object, vectors, n)

    return(object)
})
