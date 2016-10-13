#' The Params virtual class
#'
#' Virtual S4 class that all other Params classes inherit from.
#'
#' @section Parameters:
#'
#' The Params class defines the following parameters:
#'
#' \describe{
#'   \item{\code{[nGenes]}}{The number of genes to simulate.}
#'   \item{\code{[nCells]}}{The number of cells to simulate.}
#'   \item{\code{seed}}{Seed to use for generating random numbers.}
#' }
#'
#' The parameters shown in brackets can be estimated from real data.
#'
#' @name Params
#' @rdname Params
#' @aliases Params-class
setClass("Params",
         contains = "VIRTUAL",
         slots = c(nGenes = "numeric",
                   nCells = "numeric",
                   seed = "numeric"),
         prototype = prototype(nGenes = 10000, nCells = 100,
                               seed = sample(1:1e6, 1)))

#' The SimpleParams class
#'
#' S4 class that holds parameters for the simple simulation.
#'
#' @section Parameters:
#'
#' The simple simulation uses the following parameters:
#'
#' \describe{
#'   \item{\code{[nGenes]}}{The number of genes to simulate.}
#'   \item{\code{[nCells]}}{The number of cells to simulate.}
#'   \item{\code{seed}}{Seed to use for generating random numbers.}
#'   \item{\code{[mean.shape]}}{The shape parameter for the mean gamma
#'   distribution.}
#'   \item{\code{[mean.rate]}}{The rate parameter for the mean gamma
#'   distribution.}
#'   \item{\code{count.disp}}{The dispersion parameter for the counts negative
#'   binomial distribution.}
#' }
#'
#' The parameters shown in brackets can be estimated from real data using
#' \code{\link{estimateSimpleParams}}. For details of the simple simulation
#' see \code{\link{simpleSimulate}}.
#'
#' @name SimpleParams
#' @rdname SimpleParams
#' @aliases SimpleParams-class
#' @exportClass SimpleParams
setClass("SimpleParams",
         contains = "Params",
         slots = c(mean.shape = "numeric",
                   mean.rate = "numeric",
                   count.disp = "numeric"),
         prototype = prototype(mean.shape = 0.4, mean.rate = 0.3,
                               count.disp = 0.1))