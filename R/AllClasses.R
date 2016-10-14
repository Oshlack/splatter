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
#'   \item{\code{nGenes}}{The number of genes to simulate.}
#'   \item{\code{nCells}}{The number of cells to simulate.}
#'   \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'   \item{\code{mean.shape}}{The shape parameter for the mean gamma
#'   distribution.}
#'   \item{\code{mean.rate}}{The rate parameter for the mean gamma
#'   distribution.}
#'   \item{\code{[count.disp]}}{The dispersion parameter for the counts negative
#'   binomial distribution.}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{simpleEstimate}}. For details of the simple simulation
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

#' The SplatParams class
#'
#' S4 class that holds parameters for the Splatter simulation.
#'
#' @section Parameters:
#'
#' The Splatter simulation requires the following parameters:
#'
#' \describe{
#'   \item{\code{nGenes}}{The number of genes to simulate.}
#'   \item{\code{nCells}}{The number of cells to simulate.}
#'   \item{\code{[groupCells]}}{Vector giving the number of cells in each
#'   simulation group/path.}
#'   \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'   \item{\emph{Mean parameters}}{
#'     \describe{
#'       \item{\code{mean.shape}}{Shape parameter for the mean gamma
#'       distribution.}
#'       \item{\code{mean.rate}}{Rate parameter for the mean gamma
#'       distribution.}
#'     }
#'   }
#'   \item{\emph{Library size parameters}}{
#'     \describe{
#'       \item{\code{lib.loc}}{Location (meanlog) parameter for the library
#'       size log-normal distribution.}
#'       \item{\code{lib.scale}}{Scale (sdlog) parameter for the library size
#'       log-normal distribution.}
#'     }
#'   }
#'   \item{\emph{Expression outlier parameters}}{
#'     \describe{
#'       \item{\code{out.prob}}{Probability that a gene is an expression
#'       outlier.}
#'       \item{\code{out.loProb}}{Probability that an expression outlier gene
#'       is lowly expressed.}
#'       \item{\code{out.facLoc}}{Location (meanlog) parameter for the
#'       expression outlier factor log-normal distribution.}
#'       \item{\code{out.facScale}}{Scale (sdlog) parameter for the expression
#'       outlier factor log-normal distribution.}
#'     }
#'   }
#'   \item{\emph{Differential expression parameters}}{
#'     \describe{
#'       \item{\code{[de.prob]}}{Probability that a gene is differentially
#'       expressed in a group. Can be a vector.}
#'       \item{\code{[de.loProb]}}{Probability that differentially expressed
#'       gene is down-regulated. Can be a vector.}
#'       \item{\code{[de.facLoc]}}{Location (meanlog) parameter for the
#'       differential expression factor log-normal distribution. Can be a
#'       vector.}
#'       \item{\code{[de.facScale]}}{Scale (sdlog) parameter for the
#'       differential expression factor log-normal distribution. Can be a
#'       vector.}
#'     }
#'   }
#'   \item{\emph{Biological Coefficient of Variation parameters}}{
#'     \describe{
#'       \item{\code{bcv.common}}{Underlying common dispersion across all
#'       genes.}
#'       \item{\code{bcv.df}}{Degrees of Freedom for the BCV inverse chi-squared
#'       distribution.}
#'     }
#'   }
#'   \item{\emph{Dropout parameters}}{
#'     \describe{
#'       \item{\code{dropout.present}}{Logical. Whether to simulate dropout.}
#'       \item{\code{dropout.mid}}{Midpoint parameter for the dropout logistic
#'       function.}
#'       \item{\code{dropout.shape}}{Shape parameter for the dropout logistic
#'       function.}
#'     }
#'   }
#'   \item{\emph{Differentiation path parameters}}{
#'     \describe{
#'       \item{\code{[path.from]}}{Vector giving the originating point of each
#'       path. This allows path structure such as a cell type which
#'       differentiates into an intermediate cell type that then differentiates
#'       into two mature cell types. A path structure of this form would have a
#'       "from" parameter of c(0, 1, 1) (where 0 is the origin). If no vector is
#'       given all paths will start at the origin.}
#'       \item{\code{[path.length]}}{Vector giving the number of steps to
#'       simulate along each path. If a single value is given it will be applied
#'       to all paths.}
#'       \item{\code{[path.skew]}}{Vector giving the skew of each path. Values
#'       closer to 1 will give more cells towards the starting population,
#'       values closer to 0 will give more cells towards the final population.
#'       If a single value is given it will be applied to all paths.}
#'       \item{\code{[path.nonlinearProb]}}{Probability that a gene follows a
#'       non-linear path along the differentiation path. This allows more
#'       complex gene patterns such as a gene being equally expressed at the
#'       beginning an end of a path but lowly expressed in the middle.}
#'       \item{\code{[path.sigmaFac]}}{Sigma factor for non-linear gene paths.
#'       A higher value will result in more extreme non-linear variations along
#'       a path.}
#'     }
#'   }
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{splatEstimate}}. For details of the Splatter simulation
#' see \code{\link{splatSimulate}}.
#'
#' @name SplatParams
#' @rdname SplatParams
#' @aliases SplatParams-class
#' @exportClass SplatParams
setClass("SplatParams",
         contains = "Params",
         slots = c(nGroups = "numeric",
                   groupCells = "numeric",
                   mean.shape = "numeric",
                   mean.rate = "numeric",
                   lib.loc = "numeric",
                   lib.scale = "numeric",
                   out.prob = "numeric",
                   out.loProb = "numeric",
                   out.facLoc = "numeric",
                   out.facScale = "numeric",
                   de.prob = "numeric",
                   de.downProb = "numeric",
                   de.facLoc = "numeric",
                   de.facScale = "numeric",
                   bcv.common = "numeric",
                   bcv.df = "numeric",
                   dropout.present = "logical",
                   dropout.mid = "numeric",
                   dropout.shape = "numeric",
                   path.from = "numeric",
                   path.length = "numeric",
                   path.skew = "numeric",
                   path.nonlinearProb = "numeric",
                   path.sigmaFac = "numeric"),
         prototype = prototype(nGroups = 1,
                               groupCells = 100,
                               mean.rate = 0.3,
                               mean.shape = 0.4,
                               lib.loc = 10,
                               lib.scale = 0.5,
                               out.prob = 0.1,
                               out.loProb = 0.5,
                               out.facLoc = 4,
                               out.facScale = 1,
                               de.prob = 0.1,
                               de.downProb = 0.5,
                               de.facLoc = 4,
                               de.facScale = 1,
                               bcv.common = 0.1,
                               bcv.df = 25,
                               dropout.present = TRUE,
                               dropout.mid = 0,
                               dropout.shape = -1,
                               path.from = 0,
                               path.length = 100,
                               path.skew = 0.5,
                               path.nonlinearProb = 0.1,
                               path.sigmaFac = 0.8))