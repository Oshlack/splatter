#' New Params
#'
#' Create a new Params object. Functions exist for each of the different
#' Params subtypes.
#'
#' @param ... additional parameters passed to \code{\link{setParams}}.
#'
#' @return New Params object.
#' @examples
#' params <- newSimpleParams()
#' params <- newSimpleParams(nGenes = 200, nCells = 10)
#'
#' @name newParams
NULL

#' Get a parameter
#'
#' Accessor function for getting parameter values.
#'
#' @param object object to get parameter from.
#' @param name name of the parameter to get.
#'
#' @return The extracted parameter value
#'
#' @examples
#' params <- newSimpleParams()
#' getParam(params, "nGenes")
#'
#' @rdname getParam
#' @export
setGeneric("getParam", function(object, name) {standardGeneric("getParam")})

#' Set a parameter
#'
#' Function for setting parameter values.
#'
#' @param object object to set parameter in.
#' @param name name of the parameter to set.
#' @param value value to set the parameter to.
#'
#' @return Object with new parameter value.
#'
#' @examples
#' params <- newSimpleParams()
#' setParam(params, "nGenes", 100)
#'
#' @rdname setParam
#' @export
setGeneric("setParam", function(object, name, value) {
    standardGeneric("setParam")
})

#' Set a parameter UNCHECKED
#'
#' Function for setting parameter values. THE OUTPUT IS NOT CHECKED FOR
#' VALIDITY!
#'
#' @param object object to set parameter in.
#' @param name name of the parameter to set.
#' @param value value to set the parameter to.
#'
#' @return Object with new parameter value.
#'
#' @rdname setParamUnchecked
setGeneric("setParamUnchecked", function(object, name, value) {
    standardGeneric("setParamUnchecked")
})

#' Set parameters
#'
#' Set multiple parameters in a Params object.
#'
#' @param object Params object to set parameters in.
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
#'
#' @examples
#' params <- newSimpleParams()
#' params
#' # Set individually
#' params <- setParams(params, nGenes = 1000, nCells = 50)
#' params
#' # Set via update list
#' params <- setParams(params, list(mean.rate = 0.2, mean.shape = 0.8))
#' params
#'
#' @rdname setParams
#' @export
setGeneric("setParams", function(object, update = NULL, ...) {
    standardGeneric("setParams")
})

#' Expand parameters
#'
#' Expand the parameters that can be vectors so that they are the same length as
#' the number of groups.
#'
#' @param object object to expand.
#' @param ... additional arguments.
#'
#' @return Expanded object.
#'
#' @rdname expandParams
setGeneric("expandParams", function(object, ...) {
    standardGeneric("expandParams")
})
