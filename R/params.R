#' Create splatParams object
#'
#' S3 class for holding Splatter simulation parameters.
#'
#' @param ... parameters to set in the new params object, passed to
#'            \code{\link{updateParams}}.
#'
#' @details
#' The splatParams object is a list based S3 object for holding simulation
#' parameters. It has the following sections and values:
#'
#' \itemize{
#'   \item nGenes - Number of genes to simulate.
#'   \item nCells - Number of cells to simulate.
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
    params <- list(nGenes = NA, nCells = NA, groupCells = NA,
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

    params <- updateParams(params, ...)

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
#' @examples
#' params <- defaultParams()
#' print(params)
#' @export
print.splatParams <- function(x, ...) {

    pp <- list("Global:" = c("(Genes)" = x$nGenes, "(Cells)" = x$nCells,
                             "[Group Cells]" = x$groupCells),
               "Mean:" = c("(Rate)" = x$mean$rate, "(Shape)" = x$mean$shape),
               "Library size:" = c("(Location)" = x$lib$loc,
                                   "(Scale)" = x$lib$scale),
               "Expression outliers:" = c("(Probability)" = x$out$prob,
                                          "(Lo Probability)" = x$out$loProb,
                                          "(Location)" = x$out$facLoc,
                                          "(Scale)" = x$out$facScale),
               "Differential expression:" = c("[Probability]" = x$de$prob,
                                              "[Down Prob]" = x$de$downProb,
                                              "[Location]" = x$de$facLoc,
                                              "[Scale]" = x$de$facScale),
               "BCV:" = c("(Common Disp)" = x$bcv$common,
                          "(Degrees of Freedom)" = x$bcv$DF),
               "Dropout:" = c("(Present T/F)" = x$dropout$present,
                              "(Midpoint)" = x$dropout$mid,
                              "(Shape)" = x$dropout$shape),
               "Paths:" = c("[From]" = x$path$from, "[Length]" = x$path$length,
                            "[Skew]" = x$path$skew,
                            "[Non-linear Prob]" = x$path$nonlinearProb,
                            "[Sigma Factor]" = x$path$sigmaFac))

    for (category in names(pp)) {
        cat(category, "\n")
        print.default(pp[[category]], print.gap = 2)
        cat("\n")
    }
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
    types <- c(nGenes = "INT", nCells = "INT", groupCells = "INT",
               mean.rate = "POS", mean.shape = "POS", lib.loc = "NUM",
               lib.scale = "POS", out.prob = "PROB", out.loProb = "PROB",
               out.facLoc = "NUM", out.facScale = "POS", de.prob = "PROB",
               de.downProb = "PROB", de.facLoc = "NUM", de.facScale = "POS",
               bcv.common = "POS", bcv.DF = "POS", dropout.present = "LOG",
               dropout.mid = "NUM", dropout.shape = "NUM", path.from = "INT",
               path.length = "INT", path.skew = "PROB",
               path.nonlinearProb = "PROB", path.sigmaFac = "POS")

    # Define which parameters are allowed to be vectors
    vectors <- c("groupCells", "path.from", "path.length", "path.skew")
    n.groups <- length(params$groupCells)

    for (idx in seq_along(types)) {
        name <- names(types)[idx]
        name.split <- strsplit(name, ".", fixed = TRUE)[[1]]
        type <- types[idx]
        if (length(name.split) == 1) {
            value <- params[[name]]
        } else {
            value <- params[[name.split[1]]][[name.split[2]]]
        }

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

        # Missing values are allowed so we skip anything that is not NA
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
}

#' Update a splatParams object
#'
#' Update any of the parameters in a splatParams object to have a new value.
#'
#' @param params the splatParams object to update.
#' @param ... Any parameters to update.
#'
#' @details
#' This function allows multiple parameters to be updated or set using a single
#' simple function call. Parameters to update are specified by supplying
#' additional arguments that follow the levels of the splatParams data structure
#' separated by the "." character. For example
#' \code{updateParams(params, nGenes = 100)} is equivalent to
#' \code{params$nGenes <- 100} and \code{update(params, mean.rate = 1)} is
#' equivalent to \code{params$mean$rate <- 1}. For more details of the available
#' parameters and the splatParams data structure see \code{\link{splatParams}}.
#'
#' @return splatParms object with updated parameters
#' @examples
#' params <- defaultParams()
#' params
#' # Set nGenes and nCells
#' params <- updateParams(params, nGenes = 1000, nCells = 200)
#' params
#' # Set mean rate paramater and library size location parameter
#' params <- updateParams(params, mean.rate = 1, lib.loc = 12)
#' params
#' @export
updateParams <- function(params, ...) {

    update <- list(...)

    if (length(update) == 0) {
        return(params)
    }

    update.names <- strsplit(names(update), ".", fixed = TRUE)

    for (idx in 1:length(update)) {
        update.name <- update.names[[idx]]
        value <- update[[idx]]
        if (length(update.name) == 1) {
            params[[update.name]] <- value
        } else {
            params[[update.name[1]]][[update.name[2]]] <- value
        }
    }

    return(params)
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
            if (is.na(params1[[i]][[j]])) {
                params1[[i]][[j]] <- params2[[i]][[j]]
            }
        }
    }

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

    params <- updateParams(params, nGenes = 10000, nCells = 100,
                           groupCells = 100, mean.rate = 0.3, mean.shape = 0.4,
                           lib.loc = 10, lib.scale = 0.5, out.prob = 0.1,
                           out.loProb = 0.5, out.facLoc = 4, out.facScale = 1,
                           de.prob = 0.1, de.downProb = 0.5, de.facLoc = 4,
                           de.facScale = 1, bcv.common = 0.1, bcv.DF = 25,
                           dropout.present = TRUE, dropout.mid = 0,
                           dropout.shape = -1, path.from = 0,
                           path.length = 100, path.skew = 0.5,
                           path.nonlinearProb = 0.1, path.sigmaFac = 0.8)

    return(params)
}