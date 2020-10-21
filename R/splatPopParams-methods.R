#' @rdname newParams
#' @importFrom methods new
#' @export
newSplatPopParams <- function(...) {

    params <- new("SplatPopParams")
    params <- setParams(params, ...)

    checkmate::assertClass(params, "SplatPopParams")

    for (pkg in c("VariantAnnotation", "preprocessCore")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("The splatPop simulation requires the ", pkg, " package.")
        }
    }

    return(params)
}


#' @importFrom checkmate checkInt checkIntegerish checkNumber checkNumeric
#' checkFlag
setValidity("SplatPopParams", function(object) {

    v <- getParams(object, c(slotNames(object)))

    checks <- c(eqtl.n = checkNumber(v$eqtl.n, lower = 0),
                eqtl.dist = checkInt(v$eqtl.dist, lower = 1),
                eqtl.maf.min = checkNumber(v$eqtl.maf.min, lower = 0,
                                           upper = 0.5),
                eqtl.maf.max = checkNumber(v$eqtl.maf.max, lower = 0,
                                           upper = 0.5),
                eqtl.ES.shape = checkNumber(v$eqtl.ES.shape, lower = 0),
                eqtl.ES.rate = checkNumber(v$eqtl.ES.rate, lower = 0),
                eqtl.group.specific = checkNumber(v$eqtl.group.specific,
                                                 lower = 0, upper = 1),
                pop.mean.shape = checkNumber(v$pop.mean.shape, lower = 0),
                pop.mean.rate = checkNumber(v$pop.mean.rate, lower = 0),
                pop.cv.bins = checkInt(v$pop.cv.bins, lower = 1),
                pop.cv.param = checkDataFrame(v$pop.cv.param),
                similarity.scale = checkNumber(v$similarity.scale, lower = 0))

    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})


#' @importFrom methods callNextMethod
setMethod("show", "SplatPopParams", function(object) {

    pp <- list("Population params:" = c("(mean.shape)" = "pop.mean.shape",
                                        "(mean.rate)" = "pop.mean.rate",
                                        #"[group.prop]" = "group.prop",
                                        "[similarity.scale]" = "similarity.scale",
                                        "[cv.bins]" = "pop.cv.bins",
                                        "(cv.params)" = "pop.cv.param"),
               "eQTL params:" = c("[eqtl.n]"    = "eqtl.n",
                                  "[eqtl.dist]" = "eqtl.dist",
                                  "[eqtl.maf.min]" = "eqtl.maf.min",
                                  "[eqtl.maf.max]" = "eqtl.maf.max",
                                  "[eqtl.group.specific]" = "eqtl.group.specific",
                                  "(eqtl.ES.shape)" = "eqtl.ES.shape",
                                  "(eqtl.ES.rate)" = "eqtl.ES.rate"))

    callNextMethod()
    showPP(object, pp)
})


#' @rdname setParam
setMethod("setParam", "SplatPopParams", function(object, name, value) {
    checkmate::assertString(name)

    # splatPopParam checks
    if (name == "pop.cv.param") {
        if (getParam(object, "pop.cv.bins") != nrow(value)) {
            stop("Need to set pop.cv.bins to length of pop.cv.param")
        }
    }

    if (name == "eqtl.maf.min") {
        if (getParam(object, "eqtl.maf.min") >= getParam(object, "eqtl.maf.max")) {
            stop("Range of acceptable Minor Allele Frequencies is too small...
                 Be sure eqtl.maf.min < eqtl.maf.max.")
        }
    }

    # splatParam checks
    if (name == "path.length") {
        warning("path.length has been renamed path.nSteps, ",
                "please use path.nSteps in the future.")
        name <- "path.nSteps"
    }

    if (name == "nCells" || name == "nBatches") {
        stop(name, " cannot be set directly, set batchCells instead")
    }

    if (name == "nGroups") {
        stop(name, " cannot be set directly, set group.prob instead")
    }

    if (name == "batchCells") {
        object <- setParamUnchecked(object, "nCells", sum(value))
        object <- setParamUnchecked(object, "nBatches", length(value))
    }

    if (name == "group.prob") {
        object <- setParamUnchecked(object, "nGroups", length(value))
        path.from <- getParam(object, "path.from")
        if (length(path.from) > 1 & length(path.from) != length(value)) {
            warning("nGroups has changed, resetting path.from")
            object <- setParam(object, "path.from", 0)
        }
    }

    if (name == "dropout.type") {
        mid.len <- length(getParam(object, "dropout.mid"))
        mid.shape <- length(getParam(object, "dropout.shape"))
        if ((value == "experiment")) {
            if ((mid.len != 1) | (mid.shape != 1)) {
                stop("dropout.type cannot be set to 'experiment' because ",
                     "dropout.mid and dropout.shape aren't length 1, ",
                     "set dropout.mid and dropout.shape first")
            }
        }
        if ((value == "batch")) {
            n <- getParam(object, "nBatches")
            if ((mid.len != n) | (mid.shape != n)) {
                stop("dropout.type cannot be set to 'batch' because ",
                     "dropout.mid and dropout.shape aren't length equal to ",
                     "nBatches (", n, "), set dropout.mid and dropout.shape ",
                     "first")
            }
        }
        if ((value == "group")) {
            n <- getParam(object, "nGroups")
            if ((mid.len != n) | (mid.shape != n)) {
                stop("dropout.type cannot be set to 'group' because ",
                     "dropout.mid and dropout.shape aren't length equal to ",
                     "nGroups (", n, "), set dropout.mid and dropout.shape ",
                     "first")
            }
        }
        if ((value == "cell")) {
            n <- getParam(object, "nCells")
            if ((mid.len != n) | (mid.shape != n)) {
                stop("dropout.type cannot be set to 'cell' because ",
                     "dropout.mid and dropout.shape aren't length equal to ",
                     "nCells (", n, "), set dropout.mid and dropout.shape ",
                     "first")
            }
        }
    }

    object <- callNextMethod()

    return(object)
})

