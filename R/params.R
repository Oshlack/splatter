#' Create splatParams object
#'
#' S3 class for holding Splatter simulation parameters.
#'
#' @param ... parameters to set in the new params object, passed to
#'            \code{\link{setParams}}.
#'
#' @details
#' The splatParams object is a list based S3 object for holding simulation
#' parameters. It has the following sections and values:
#'
#' \itemize{
#'   \item nGenes - Number of genes to simulate.
#'   \item nCells - Number of cells to simulate.
#'   \item nGroups - Number of groups to simulate.
#'   \item [groupCells] - Vector giving the number of cells in each simulation
#'         group/path.
#'   \item mean (mean parameters)
#'     \itemize{
#'       \item rate - Rate parameter for the mean gamma simulation.
#'       \item shape - Shape parameter for the mean gamma simulation.
#'     }
#'   \item lib (library size parameters)
#'     \itemize{
#'       \item loc - Location (meanlog) parameter for the library size
#'             log-normal distribution.
#'       \item scale - Scale (sdlog) parameter for the library size log-normal
#'             distribution.
#'     }
#'   \item out (expression outlier parameters)
#'     \itemize{
#'       \item prob - Probability that a gene is an expression outlier.
#'       \item loProb - Probability that an expression outlier gene is lowly
#'             expressed.
#'       \item facLoc - Location (meanlog) parameter for the expression outlier
#'             factor log-normal distribution.
#'       \item facScale - Scale (sdlog) parameter for the expression outlier
#'             factor log-normal distribution.
#'     }
#'   \item de (differential expression parameters)
#'     \itemize{
#'       \item [prob] - Probability that a gene is differentially expressed
#'             between groups or paths.
#'       \item [downProb] - Probability that differentially expressed gene is
#'             down-regulated.
#'       \item [facLoc] - Location (meanlog) parameter for the differential
#'             expression factor log-normal distribution.
#'       \item [facScale] - Scale (sdlog) parameter for the differential
#'             expression factor log-normal distribution.
#'     }
#'   \item bcv (Biological Coefficient of Variation parameters)
#'     \itemize{
#'       \item common - Underlying common dispersion across all genes.
#'       \item DF - Degrees of Freedom for the BCV inverse chi-squared
#'             distribution.
#'     }
#'   \item dropout (dropout parameters)
#'     \itemize{
#'       \item present - Logical. Whether to simulate dropout.
#'       \item mid - Midpoint parameter for the dropout logistic function.
#'       \item shape - Shape parameter for the dropout logistic function.
#'     }
#'   \item path (differentiation path parameters)
#'     \itemize{
#'       \item [from] - Vector giving the originating point of each path. This
#'             allows path structure such as a cell type which differentiates
#'             into an intermediate cell type that then differentiates into two
#'             mature cell types. A path structure of this form would have a
#'             "from" parameter of c(0, 1, 1) (where 0 is the origin). If no
#'             vector is given all paths will start at the origin.
#'       \item [length] - Vector giving the number of steps to simulate along
#'             each path. If a single value is given it will be applied to all
#'             paths.
#'       \item [skew] - Vector giving the skew of each path. Values closer to 1
#'             will give more cells towards the starting population, values
#'             closer to 0 will give more cells towards the final population.
#'             If a single value is given it will be applied to all paths.
#'       \item [nonlinearProb] - Probability that a gene follows a non-linear
#'             path along the differentiation path. This allows more complex
#'             gene patterns such as a gene being equally expressed at the
#'             beginning an end of a path but lowly expressed in the middle.
#'       \item [sigmaFac] - Sigma factor for non-linear gene paths. A higher
#'             value will result in more extreme non-linear variations along a
#'             path.
#'   }
#' }
#'
#' Those shown in brackets cannot currently be easily estimated from a real
#' dataset by Splatter. This is also shown when a splatParams object is printed
#' with parameters labelled as either (estimatable) or [not estimatable].
#'
#' @return List based S3 splatParams object
#' @examples
#' params <- splatParams()
#' params
#' @export
splatParams <- function(...) {
    params <- list(nGenes = NA, nCells = NA, nGroups = NA, groupCells = NA,
                   mean = list(rate = NA, shape = NA),
                   lib = list(loc = NA, scale = NA),
                   out = list(prob = NA, loProb = NA, facLoc = NA,
                              facScale = NA),
                   de = list(prob = NA, downProb = NA, facLoc = NA,
                             facScale = NA),
                   bcv = list(common = NA, DF = NA),
                   dropout = list(present = NA, mid = NA, shape = NA),
                   path = list(from = NA, length = NA, skew = NA,
                               nonlinearProb = NA, sigmaFac = NA))

    class(params) <- "splatParams"

    params <- setParams(params, ...)

    return(params)
}

#' Print splatParams object
#'
#' Pretty print the parameters stored in a splatParams object. Parameters are
#' labelled as either (estimatable) or [not estimatable].
#'
#' @param x splatParams object to print.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Invisibly returns x (via \code{\link{invisible}(x)})
#'
#' @examples
#' params <- defaultParams()
#' print(params)
#' @export
print.splatParams <- function(x, ...) {

    pp <- list("Global:"             = c("(Genes)"        = "nGenes",
                                         "(Cells)"        = "nCells",
                                         "[Groups]"       = "nGroups",
                                         "[Group Cells]"  = "groupCells"),
               "Mean:"               = c("(Rate)"         = "mean.rate",
                                         "(Shape)"        = "mean.shape"),
               "Library size:"       = c("(Location)"     = "lib.loc",
                                         "(Scale)"        = "lib.scale"),
               "Exprs outliers:"     = c("(Probability)"  = "out.prob",
                                         "(Lo Prob)"      = "out.loProb",
                                         "(Location)"     = "out.facLoc",
                                         "(Scale)"        = "out.facScale"),
               "Differential exprs:" = c("[Probability]"  = "de.prob",
                                         "[Down Prob]"    = "de.downProb",
                                         "[Location]"     = "de.facLoc",
                                         "[Scale]"        = "de.facScale"),
               "BCV:"                = c("(Common Disp)"  = "bcv.common",
                                         "(DoF)"          = "bcv.DF"),
               "Dropout:"            = c("(Present T/F)"  = "dropout.present",
                                         "(Midpoint)"     = "dropout.mid",
                                         "(Shape)"        = "dropout.shape"),
               "Paths:"              = c("[From]"         = "path.from",
                                         "[Length]"       = "path.length",
                                         "[Skew]"         = "path.skew",
                                         "[Non-linear]"   = "path.nonlinearProb",
                                         "[Sigma Factor]" = "path.sigmaFac"))

    for (category in names(pp)) {
        parameters <- getParams(x, pp[[category]])
        parameters <- sapply(parameters, paste, collapse = ", ")
        names(parameters) <- names(pp[[category]])
        cat(category, "\n")
        print(noquote(parameters), print.gap = 2)
        cat("\n")
    }

    invisible(x)
}

#' Update a splatParams object
#'
#' Set any of the parameters in a splatParams object to have a new value.
#'
#' @param params the splatParams object to update.
#' @param ... Any parameters to set.
#'
#' @details
#' This function allows multiple parameters to be updated or set using a single
#' simple function call. Parameters to update are specified by supplying
#' additional arguments that follow the levels of the splatParams data structure
#' separated by the "." character. For example
#' \code{setParams(params, nGenes = 100)} is equivalent to
#' \code{params$nGenes <- 100} and \code{update(params, mean.rate = 1)} is
#' equivalent to \code{params$mean$rate <- 1}. For more details of the available
#' parameters and the splatParams data structure see \code{\link{splatParams}}.
#'
#' @return splatParms object with updated parameters
#' @examples
#' params <- defaultParams()
#' params
#' # Set nGenes and nCells
#' params <- setParams(params, nGenes = 1000, nCells = 200)
#' params
#' # Set mean rate paramater and library size location parameter
#' params <- setParams(params, mean.rate = 1, lib.loc = 12)
#' params
#' @export
setParams <- function(params, ...) {

    update <- list(...)

    if (length(update) == 0) {
        return(params)
    }

    names <- names(update)

    for (idx in seq_along(names)) {
        name <- names[idx]
        name.split <- strsplit(name, ".", fixed = TRUE)[[1]]
        value <- update[[idx]]

        if (name == "nCells" || name == "nGroups") {
            stop(name, " cannot be set directly, set groupCells instead")
        }

        if (length(name.split) == 1) {
            params[[name.split]] <- value
        } else {
            params[[name.split[1]]][[name.split[2]]] <- value
        }

        if (name == "groupCells") {
            params$nCells <- sum(value)
            params$nGroups <- length(value)
        }
    }

    checkParams(params)

    return(params)
}

#' Get parameters from splatParams object
#'
#' Get values for the parameters in a splatParams object. Uses the same pattern
#' (category.parameter) as \code{\link{setParams}}.
#'
#' @param params splatParams object to get parameters from.
#' @param names vector of parameter names to extract.
#'
#' @return Vector if all selected parameters are single values, otherwise a
#'         list.
#' @examples
#' params <- defaultParams()
#' # Get the number of genes
#' getParams(params, "nGenes")
#' # Get the number of genes and the mean rate parameter
#' getParams(params, c("nGenes", "mean.rate"))
#' # Returns a list if one of the selected parameters is a vector
#' params <- setParams(params, groupCells = c(100, 200))
#' getParams(params, c("nGenes", "mean.rate", "groupCells"))
#' @export
getParams <- function(params, names) {

    if (length(names) == 0) {
        return(NULL)
    }

    output <- list()
    keep.list <- FALSE
    for (idx in seq_along(names)) {
        name <- names[idx]
        name.split <- strsplit(name, ".", fixed = TRUE)[[1]]
        if (length(name.split) == 1) {
            value <- params[[name.split]]
        } else {
            value <- params[[name.split[1]]][[name.split[2]]]
        }
        output[[name]] <- value
        if (length(value) > 1) {
            keep.list <- TRUE
        }
    }

    if (!keep.list) {
        output <- unlist(output)
    } else if (length(output) == 1) {
        output <- output[[1]]
    }

    return(output)
}

#' Check splatParams object
#'
#' Check that a splatParams object has valid parameter values.
#'
#' @param params splatParams object to check
#'
#' @details
#' The following checks are made:
#' \itemize{
#'   \item{Input has "splatParams" class}
#'   \item{Logical parameters are logical}
#'   \item{Numeric parameters are numeric}
#'   \item{Positive numeric parameters are positive}
#'   \item{Integer parameters are integers}
#'   \item{Probability parameters are in the range 0-1}
#'   \item{Vector parameters are the correct length}
#'   \item{Vector parameters do not contain NAs}
#'   \item{Non-vector parameters are single values}
#'   \item{nCells, nGroups and groupCells are consistent}
#'   \item{path.from has possible values}
#' }
#'
#' @return Produces error if not valid otherwise nothing
#' @examples
#' checkParams(defaultParams())
#' @export
checkParams <- function(params) {

    # Check class before anything else
    if (!("splatParams" %in% class(params))) {
        stop("params does not belong to the splatParams class")
    }

    # Define what values each parameter can take
    # NUM = Numeric
    # POS = Positive numeric
    # INT = Positive integer
    # PROB = Positive numeric in range 0-1
    # LOG = Logical
    types <- c(nGenes = "INT", nCells = "INT", nGroups = "INT",
               groupCells = "INT", mean.rate = "POS", mean.shape = "POS",
               lib.loc = "NUM", lib.scale = "POS", out.prob = "PROB",
               out.loProb = "PROB", out.facLoc = "NUM", out.facScale = "POS",
               de.prob = "PROB", de.downProb = "PROB", de.facLoc = "NUM",
               de.facScale = "POS", bcv.common = "POS", bcv.DF = "POS",
               dropout.present = "LOG", dropout.mid = "NUM",
               dropout.shape = "NUM", path.from = "INT", path.length = "INT",
               path.skew = "PROB", path.nonlinearProb = "PROB",
               path.sigmaFac = "POS")

    # Define which parameters are allowed to be vectors
    vectors <- c("groupCells", "path.from", "path.length", "path.skew")
    n.groups <- length(getParams(params, "groupCells"))

    for (idx in seq_along(types)) {
        name <- names(types)[idx]
        type <- types[idx]
        value <- getParams(params, name)

        # Check vector properties first so we can exclude vectors with an NA
        # before the next section
        if (length(value) > 1) {
            if (name %in% vectors) {
                if (any(is.na(value))) {
                    stop(name, " is a vector and contains NA values")
                } else if (length(value) != n.groups) {
                    stop("length of ", name, " must be 1 or the length of ",
                         "the groupCells parameter")
                }
            } else {
                stop(name, " should be a single value")
            }
        }

        # Missing values are allowed so we skip anything that is NA
        if (!all(is.na(value))) {

            if (type %in% c("NUM", "INT", "POS", "PROB") &&
                !(is.numeric(value))) {
                stop(name, " must be numeric")
            }

            if (type %in% c("INT", "POS", "PROB") && value < 0) {
                stop(name, " must be positive")
            }

            if (type == "INT" && value %% 1 != 0) {
                stop(name, " must be an integer")
            }

            if (type == "PROB" && (value < 0 || value > 1)) {
                stop(paste(name, "must be in the range 0-1"))
            }

            if (type == "LOG" && !(is.logical(value))) {
                stop(name, " must be logical (TRUE/FALSE)")
            }
        }
    }

    # Check groupCells matches nCells, nGroups
    n.cells <- getParams(params, "nCells")
    n.groups <- getParams(params, "nGroups")
    group.cells <- getParams(params, "groupCells")
    if (!all(is.na(group.cells)) &&
        (n.cells != sum(group.cells) || n.groups != length(group.cells))) {
        stop("nCells, nGroups and groupCells are not consistent")
    }

    # Check path.from contains origin, has allowed values and does not reference
    # itself (no loops)
    path.from <- getParams(params, "path.from")
    if (!all(is.na(path.from))) {
        if (!(0 %in% path.from)) {
            stop("origin must be specified in path.from")
        } else if (any(path.from > n.groups)) {
            stop("values in path.from cannot be greater than number of paths")
        } else if (any(path.from == 1:n.groups)) {
            stop("path cannot begin at itself")
        }
    }
}

#' Merge two splatParams objects
#'
#' Merge two splatParams objects. Any parameters that are NA in the first
#' splatParams object are replaced by the value in the second splatParams
#' object.
#'
#' @param params1 first splatParams object to merge.
#' @param params2 second splatParams object to merge.
#'
#' @return Merged splatParams object.
#' @examples
#' params <- splatParams(nGenes = 1000, nCells = 50)
#' params
#' # Replace unset parameters with default parameters
#' params <- mergeParams(params, defaultParams())
#' params
#' @export
mergeParams <- function(params1, params2) {

    for (i in 1:length(params1)) {
        for (j in 1:length(params1[[i]])) {
            if (all(is.na(params1[[i]][[j]]))) {
                params1[[i]][[j]] <- params2[[i]][[j]]
            }
        }
    }

    checkParams(params1)

    return(params1)
}

#' Get default simulation parameters
#'
#' Get a splatParams object with a set of default parameters that will produce a
#' resonable simulation of single-cell RNA-seq count data.
#'
#' @return A splatParams object containing default parameters
#' @examples
#' params <- defaultParams()
#' params
#' @export
defaultParams <- function() {

    params <- splatParams()

    params <- setParams(params, nGenes = 10000, groupCells = 100,
                        mean.rate = 0.3, mean.shape = 0.4,
                        lib.loc = 10, lib.scale = 0.5, out.prob = 0.1,
                        out.loProb = 0.5, out.facLoc = 4, out.facScale = 1,
                        de.prob = 0.1, de.downProb = 0.5, de.facLoc = 4,
                        de.facScale = 1, bcv.common = 0.1, bcv.DF = 25,
                        dropout.present = TRUE, dropout.mid = 0,
                        dropout.shape = -1, path.from = 0, path.length = 100,
                        path.skew = 0.5, path.nonlinearProb = 0.1,
                        path.sigmaFac = 0.8)

    return(params)
}


#' Expath path parameters
#'
#' Expand the path parameters so that they are the same length as the number
#' of groups.
#'
#' @param params splatParams object to expand.
#'
#' @return expanded splatParams object.
#' @examples
#' params <- defaultParams()
#' params <- setParams(params, groupCells = c(10, 10))
#' params
#' params <- expandPathParams(params)
#' params
expandPathParams <- function(params) {

    n.groups <- getParams(params, "nGroups")
    path.from <- getParams(params, "path.from")
    path.length <- getParams(params, "path.length")
    path.skew <- getParams(params, "path.skew")

    if (length(path.from) == 1) {
        params <- setParams(params, path.from = rep(path.from, n.groups))
    }

    if (length(path.length) == 1) {
        params <- setParams(params, path.length = rep(path.length, n.groups))
    }

    if (length(path.skew) == 1) {
        params <- setParams(params, path.skew = rep(path.skew, n.groups))
    }

    return(params)
}