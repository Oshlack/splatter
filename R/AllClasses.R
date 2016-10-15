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
#'   \item{\code{[nGroups]}}{The number of groups or paths to simulate.}
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

#' The LunParams class
#'
#' S4 class that holds parameters for the Lun simulation.
#'
#' @section Parameters:
#'
#' The Lun simulation uses the following parameters:
#'
#' \describe{
#'   \item{\code{nGenes}}{The number of genes to simulate.}
#'   \item{\code{nCells}}{The number of cells to simulate.}
#'   \item{\code{[nGroups]}}{The number of groups to simulate.}
#'   \item{\code{[groupCells]}}{Vector giving the number of cells in each
#'   simulation group/path.}
#'   \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'   \item{\emph{Mean parameters}}{
#'     \describe{
#'       \item{\code{[mean.shape]}}{Shape parameter for the mean gamma
#'       distribution.}
#'       \item{\code{[mean.rate]}}{Rate parameter for the mean gamma
#'       distribution.}
#'     }
#'   }
#'   \item{\emph{Counts parameters}}{
#'     \describe{
#'       \item{\code{[count.disp]}}{The dispersion parameter for the counts
#'       negative binomial distribution.}
#'     }
#'   }
#'   \item{\emph{Differential expression parameters}}{
#'     \describe{
#'       \item{\code{[de.nGenes]}}{The number of genes that are differentially
#'       expressed in each group}
#'       \item{\code{[de.upProb]}}{The proportion of differentially expressed
#'       genes that are up-regulated in each group}
#'       \item{\code{[de.upFC]}}{The fold change for up-regulated genes}
#'       \item{\code{[de.downFC]}}{The fold change for down-regulated genes}
#'     }
#'   }
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{lunEstimate}}. For details of the Lun simulation see
#' \code{\link{lunSimulate}}.
#'
#' @name LunParams
#' @rdname LunParams
#' @aliases LunParams-class
#' @exportClass LunParams
setClass("LunParams",
         contains = "SimpleParams",
         slots = c(nGroups = "numeric",
                   groupCells = "numeric",
                   de.nGenes = "numeric",
                   de.upProp = "numeric",
                   de.upFC = "numeric",
                   de.downFC = "numeric"),
         prototype = prototype(nGroups = 1, groupCells = 100, mean.shape = 2,
                               mean.rate = 2, de.nGenes = 1000, de.upProp = 0.5,
                               de.upFC = 5, de.downFC = 0))

#' The Lun2Params class
#'
#' S4 class that holds parameters for the Lun simulation.
#'
#' @section Parameters:
#'
#' The Lun2 simulation uses the following parameters:
#'
#' \describe{
#'   \item{\code{nGenes}}{The number of genes to simulate.}
#'   \item{\code{nCells}}{The number of cells to simulate.}
#'   \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'   \item{\code{[nPlates]}}{The number of plates to simulate.}
#'   \item{\emph{Plate parameters}}{
#'     \describe{
#'       \item{\code{plate.ingroup}}{Character vecotor giving the plates
#'       considered to be part of the "ingroup".}
#'       \item{\code{plate.mod}}{Plate effect modifier factor. The plate effect
#'       variance is divided by this value.}
#'       \item{\code{plate.var}}{Plate effect variance.}
#'     }
#'   }
#'   \item{\emph{Gene parameters}}{
#'     \describe{
#'       \item{\code{gene.means}}{Mean expression for each gene.}
#'       \item{\code{gene.disps}}{Dispersion for each gene.}
#'       \item{\code{gene.ziMeans}}{Zero-inflated gene means.}
#'       \item{\code{gene.ziDisps}}{Zero-inflated gene dispersions.}
#'       \item{\code{gene.ziProps}}{Zero-inflated gene zero proportions.}
#'     }
#'   }
#'   \item{\emph{Cell parameters}}{
#'     \describe{
#'       \item{\code{cell.plates}}{Factor giving the plate that each cell comes
#'       from.}
#'       \item{\code{cell.libSizes}}{Library size for each cell.}
#'       \item{\code{cell.libMod}}{Modifier factor for library sizes.
#'       The library sizes are multiplied by this value.}
#'     }
#'   }
#'   \item{\emph{Differential expression parameters}}{
#'     \describe{
#'       \item{\code{de.nGenes}}{Number of differentially expressed genes.}
#'       \item{\code{de.fc}}{Fold change for differentially expressed genes.}
#'     }
#'   }
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{lun2Estimate}}. For details of the Lun2 simulation see
#' \code{\link{lun2Simulate}}.
#'
#' @name Lun2Params
#' @rdname Lun2Params
#' @aliases Lun2Params-class
#' @exportClass Lun2Params
setClass("Lun2Params",
         contains = "Params",
         slots = c(nPlates = "numeric",
                   plate.ingroup = "character",
                   plate.mod = "numeric",
                   plate.var = "numeric",
                   gene.means = "numeric",
                   gene.disps = "numeric",
                   gene.ziMeans = "numeric",
                   gene.ziDisps = "numeric",
                   gene.ziProps = "numeric",
                   cell.plates = "numeric",
                   cell.libSizes = "numeric",
                   cell.libMod = "numeric",
                   de.nGenes = "numeric",
                   de.fc = "numeric"),
         prototype = prototype(nPlates = 1,
                               cell.plates = factor(rep(1, 100)),
                               plate.ingroup = "1",
                               plate.mod = 1,
                               plate.var = 14,
                               gene.means = rep(3.2, 10000),
                               gene.disps = rep(0.03, 10000),
                               gene.ziMeans = rep(1.6, 10000),
                               gene.ziDisps = rep(0.1, 10000),
                               gene.ziProps = rep(2.3e-6, 10000),
                               cell.libSizes = rep(70000, 100),
                               cell.libMod = 1,
                               de.nGenes = 0,
                               de.fc = 3))

#' The SCDDParams class
#'
#' S4 class that holds parameters for the scDD simulation.
#'
#' @section Parameters:
#'
#' The SCDD simulation uses the following parameters:
#'
#' \describe{
#'   \item{\code{nGenes}}{The number of genes to simulate.}
#'   \item{\code{nCells}}{The number of cells to simulate.}
#'   \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'   \item{\code{SCdat}}{\code{\link{ExpressSet}} containing real data.}
#'   \item{\code{[nDE]}}{Number of DE genes to simulate.}
#'   \item{\code{[nDP]}}{Number of DP genes to simulate.}
#'   \item{\code{[nDM]}}{Number of DM genes to simulate.}
#'   \item{\code{[nDB]}}{Number of DB genes to simulate.}
#'   \item{\code{[nEE]}}{Number of EE genes to simulate.}
#'   \item{\code{[nEP]}}{Number of EP genes to simulate.}
#'   \item{\code{[sd.range]}}{Interval for fold change standard deviations.}
#'   \item{\code{[modeFC]}}{Values for DP, DM and DB mode fold changes.}
#'   \item{\code{[varInflation}]}{Variance inflation factors for each
#'   condition.}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{scDDEstimate}}. See \code{\link[scDD]{simulateSet}} for more
#' details of the parameters. For details of the Splatter implementation of the
#' scDD simulation see \code{\link{scDDSimulate}}.
#'
#' @name SCDDParams
#' @rdname SCDDParams
#' @aliases SCDDParams-class
#' @exportClass SCDDParams
setClass("SCDDParams",
         contains = "Params",
         slots = c(SCDat = "ExpressionSet",
                   nDE = "numeric",
                   nDP = "numeric",
                   nDM = "numeric",
                   nDB = "numeric",
                   nEE = "numeric",
                   nEP = "numeric",
                   sd.range = "numeric",
                   modeFC = "numeric",
                   varInflation = "numeric"),
         prototype = prototype(SCDat = ExpressionSet(),
                               nCells = 100,
                               nDE = 250,
                               nDP = 250,
                               nDM = 250,
                               nDB = 250,
                               nEE = 5000,
                               nEP = 4000,
                               sd.range = c(1, 3),
                               modeFC = c(2, 3, 4),
                               varInflation = c(1, 1)))