#' Get parameters
#'
#' Get multiple parameter values from a Params object.
#'
#' @param params Params object to get values from.
#' @param names vector of names of the parameters to get.
#'
#' @return List with the values of the selected parameters.
#' @examples
#' params <- newSimpleParams()
#' getParams(params, c("nGenes", "nCells", "mean.rate"))
#' @export
getParams <- function(params, names) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertCharacter(names, min.len = 1, any.missing = FALSE)

    sapply(names, getParam, object = params, simplify = FALSE)
}

#' Set parameters
#'
#' Set multiple parameters in a Params object.
#'
#' @param params Params object to set parameters in.
#' @param update list of parameters to set where \code{names(update)} are the
#'        names of the parameters to set and the items in the list are values.
#' @param ... additional parameters to set. These are combined with any
#'        parameters specified in \code{update}.
#' @param checkValid logical. Whether to check set object is valid.
#'
#' @details
#' Each parameter is set by a call to \code{\link{setParam}}. If the same
#' parameter is specified multiple times it will be set multiple times.
#' Parameters can be specified using a list via \code{update} (useful when
#' collecting parameter values in some way) or individually (useful when setting
#' them manually), see examples.
#'
#' @return Params object with updated values.
#' @examples
#' params <- newSimpleParams()
#' params
#' # Set individually
#' params <- setParams(params, nGenes = 1000, nCells = 50)
#' params
#' # Set via update list
#' params <- setParams(params, list(mean.rate = 0.2, mean.shape = 0.8))
#' params
#' @export
setParams <- function(params, update = NULL, checkValid = TRUE, ...) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertList(update, null.ok = TRUE)

    update <- c(update, list(...))

    if (length(update) > 0) {
        for (name in names(update)) {
            value <- update[[name]]
            params <- setParam(params, name, value)
        }
    }

    return(params)
}

#' Set parameters UNCHECKED
#'
#' Set multiple parameters in a Params object.
#'
#' @param params Params object to set parameters in.
#' @param update list of parameters to set where \code{names(update)} are the
#'        names of the parameters to set and the items in the list are values.
#' @param ... additional parameters to set. These are combined with any
#'        parameters specified in \code{update}.
#' @param checkValid logical. Whether to check set object is valid.
#'
#' @details
#' Each parameter is set by a call to \code{\link{setParam}}. If the same
#' parameter is specified multiple times it will be set multiple times.
#' Parameters can be specified using a list via \code{update} (useful when
#' collecting parameter values in some way) or individually (useful when setting
#' them manually), see examples. THE FINAL OBJECT IS NOT CHECKED FOR VALIDITY!
#'
#' @return Params object with updated values.
setParamsUnchecked <- function(params, update = NULL, checkValid = TRUE, ...) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertList(update, null.ok = TRUE)

    update <- c(update, list(...))

    if (length(update) > 0) {
        for (name in names(update)) {
            value <- update[[name]]
            params <- setParamUnchecked(params, name, value)
        }
    }

    return(params)
}

#' Show pretty print
#'
#' Function used for pretty printing params object.
#'
#' @param params object to show.
#' @param pp list specifying how the object should be displayed.
#'
#' @importFrom utils head
showPP <- function(params, pp) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertList(pp, types = "character", min.len = 1)

    default <- new(class(params))
    for (category in names(pp)) {
        parameters <- pp[[category]]
        values <- getParams(params, parameters)
        short.values <- sapply(values, function(x) {
            if (length(x) > 4) {
                paste0(paste(head(x, n = 4), collapse = ", "), ",...")
            } else {
                paste(x, collapse = ", ")
            }
        })
        values <- sapply(values, paste, collapse = ", ")
        default.values <- getParams(default, parameters)
        default.values <- sapply(default.values, paste, collapse = ", ")
        not.default <- values != default.values
        names(values)[not.default] <- toupper(names(values[not.default]))
        cat(category, "\n")
        print(noquote(short.values), print.gap = 2)
        cat("\n")
    }
}

# mergeParams <- function(params1, params2) {
#
#     if (class(params1) != class(params2)) {
#         stop("params1 and params2 must be of the same Params class")
#     }
#
#     default <- new(class(params1))
#
#     update <- list()
#     for (parameter in slotNames(params1)) {
#         value1 <- getParam(params1, parameter)
#         default.value <- getParam(default, parameter)
#         if (value1 == default.value) {
#             value2 <- getParam(params2, parameter)
#             update[[parameter]] <- value2
#         } else {
#             update[[parameter]] <- value1
#         }
#     }
#
#     merged <- setParams(default, update)
#
#     return(merged)
# }