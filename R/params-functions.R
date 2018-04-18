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
setParams <- function(params, update = NULL, ...) {

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
#'
#' @details
#' Each parameter is set by a call to \code{\link{setParam}}. If the same
#' parameter is specified multiple times it will be set multiple times.
#' Parameters can be specified using a list via \code{update} (useful when
#' collecting parameter values in some way) or individually (useful when setting
#' them manually), see examples. THE FINAL OBJECT IS NOT CHECKED FOR VALIDITY!
#'
#' @return Params object with updated values.
setParamsUnchecked <- function(params, update = NULL, ...) {

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
#' @return Print params object to console
showPP <- function(params, pp) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertList(pp, types = "character", min.len = 1)

    default <- new(class(params))
    for (category in names(pp)) {
        parameters <- pp[[category]]
        values <- getParams(params, parameters)
        is.df <- sapply(values, is.data.frame)

        default.values <- getParams(default, parameters)
        not.default <- sapply(seq_along(values), function(i) {
            !identical(values[i], default.values[i])
        })

        cat(category, "\n")
        if (sum(!is.df) > 0) {
            showValues(values[!is.df], not.default[!is.df])
        }
        if (sum(is.df) > 0) {
            showDFs(values[is.df], not.default[is.df])
        }
        cat("\n")
    }
}

#' Show vales
#'
#' Function used for pretty printing scale or vector parameters.
#'
#' @param values list of values to show.
#' @param not.default logical vector giving which have changed from the default.
#'
#' @importFrom utils head
showValues <- function(values, not.default) {

    short.values <- sapply(values, function(x) {
        if (length(x) > 4) {
            paste0(paste(head(x, n = 4), collapse = ", "), ",...")
        } else {
            paste(x, collapse = ", ")
        }
    })

    names(short.values)[not.default] <- toupper(names(values[not.default]))

    print(noquote(short.values), print.gap = 2)
}

#' Show data.frame
#'
#' Function used for pretty printing data.frame parameters.
#'
#' @param dfs list of data.frames to show.
#' @param not.default logical vector giving which have changed from the default.
#'
#' @importFrom utils head
showDFs <- function(dfs, not.default) {

    names(dfs)[not.default] <- toupper(names(dfs)[not.default])

    for (i in seq_along(dfs)) {
        df <- dfs[[i]]
        name <- names(dfs)[i]

        cat(paste0("\n", name, "\n"))
        cat("data.frame", paste0("(", nrow(df), " x ", ncol(df), ")"),
            "with columns:", paste(colnames(df), collapse = ", "), "\n")
        print(head(df, n = 4))
        cat("# ... with", nrow(df) - 4, "more rows\n")
    }
}
