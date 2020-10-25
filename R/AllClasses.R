#' The Params virtual class
#'
#' Virtual S4 class that all other Params classes inherit from.
#'
#' @section Parameters:
#'
#' The Params class defines the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data.
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
                               seed = sample(seq_len(1e6), 1)))


#' The SimpleParams class
#'
#' S4 class that holds parameters for the simple simulation.
#'
#' @section Parameters:
#'
#' The simple simulation uses the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\code{mean.shape}}{The shape parameter for the mean gamma
#'     distribution.}
#'     \item{\code{mean.rate}}{The rate parameter for the mean gamma
#'     distribution.}
#'     \item{\code{[count.disp]}}{The dispersion parameter for the counts
#'     negative binomial distribution.}
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
#' S4 class that holds parameters for the Splat simulation.
#'
#' @section Parameters:
#'
#' The Splat simulation requires the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\emph{Batch parameters}}{
#'         \describe{
#'             \item{\code{[nBatches]}}{The number of batches to simulate.}
#'             \item{\code{[batchCells]}}{Vector giving the number of cells in
#'             each batch.}
#'             \item{\code{[batch.facLoc]}}{Location (meanlog) parameter for the
#'             batch effect factor log-normal distribution. Can be a vector.}
#'             \item{\code{[batch.facScale]}}{Scale (sdlog) parameter for the
#'             batch effect factor log-normal distribution. Can be a vector.}
#'             \item{\code{[batch.rmEffect]}}{Logical, removes the batch effect
#'             and continues with the simulation when TRUE. This allows the 
#'             user to test batch removal algorithms without having to calculate
#'             the new expected cell means with batch removed.}
#'         }
#'     }
#'     \item{\emph{Mean parameters}}{
#'         \describe{
#'             \item{\code{mean.shape}}{Shape parameter for the mean gamma
#'             distribution.}
#'             \item{\code{mean.rate}}{Rate parameter for the mean gamma
#'             distribution.}
#'         }
#'     }
#'     \item{\emph{Library size parameters}}{
#'         \describe{
#'             \item{\code{lib.loc}}{Location (meanlog) parameter for the
#'             library size log-normal distribution, or mean parameter if a
#'             normal distribution is used.}
#'             \item{\code{lib.scale}}{Scale (sdlog) parameter for the library
#'             size log-normal distribution, or sd parameter if a normal
#'             distribution is used.}
#'             \item{\code{lib.norm}}{Logical. Whether to use a normal
#'             distribution for library sizes instead of a log-normal.}
#'         }
#'     }
#'     \item{\emph{Expression outlier parameters}}{
#'         \describe{
#'             \item{\code{out.prob}}{Probability that a gene is an expression
#'             outlier.}
#'             \item{\code{out.facLoc}}{Location (meanlog) parameter for the
#'             expression outlier factor log-normal distribution.}
#'             \item{\code{out.facScale}}{Scale (sdlog) parameter for the
#'             expression outlier factor log-normal distribution.}
#'         }
#'     }
#'     \item{\emph{Group parameters}}{
#'         \describe{
#'             \item{\code{[nGroups]}}{The number of groups or paths to
#'             simulate.}
#'             \item{\code{[group.prob]}}{Probability that a cell comes from a
#'             group.}
#'         }
#'     }
#'     \item{\emph{Differential expression parameters}}{
#'         \describe{
#'             \item{\code{[de.prob]}}{Probability that a gene is differentially
#'             expressed in a group. Can be a vector.}
#'             \item{\code{[de.downProb]}}{Probability that a differentially
#'             expressed gene is down-regulated. Can be a vector.}
#'             \item{\code{[de.facLoc]}}{Location (meanlog) parameter for the
#'             differential expression factor log-normal distribution. Can be a
#'             vector.}
#'             \item{\code{[de.facScale]}}{Scale (sdlog) parameter for the
#'             differential expression factor log-normal distribution. Can be a
#'             vector.}
#'         }
#'     }
#'     \item{\emph{Biological Coefficient of Variation parameters}}{
#'         \describe{
#'             \item{\code{bcv.common}}{Underlying common dispersion across all
#'             genes.}
#'             \item{\code{bcv.df}}{Degrees of Freedom for the BCV inverse
#'             chi-squared distribution.}
#'         }
#'     }
#'     \item{\emph{Dropout parameters}}{
#'         \describe{
#'             \item{\code{dropout.type}}{The type of dropout to simulate.
#'             "none" indicates no dropout, "experiment" is global dropout using
#'             the same parameters for every cell, "batch" uses the same
#'             parameters for every cell in each batch, "group" uses the same
#'             parameters for every cell in each groups and "cell" uses a
#'             different set of parameters for each cell.}
#'             \item{\code{dropout.mid}}{Midpoint parameter for the dropout
#'             logistic function.}
#'             \item{\code{dropout.shape}}{Shape parameter for the dropout
#'             logistic function.}
#'         }
#'     }
#'     \item{\emph{Differentiation path parameters}}{
#'         \describe{
#'             \item{\code{[path.from]}}{Vector giving the originating point of
#'             each path. This allows path structure such as a cell type which
#'             differentiates into an intermediate cell type that then
#'             differentiates into two mature cell types. A path structure of
#'             this form would have a "from" parameter of c(0, 1, 1) (where 0 is
#'             the origin). If no vector is given all paths will start at the
#'             origin.}
#'             \item{\code{[path.nSteps]}}{Vector giving the number of steps to
#'             simulate along each path. If a single value is given it will be
#'             applied to all paths. This parameter was previously called
#'             \code{path.length}.}
#'             \item{\code{[path.skew]}}{Vector giving the skew of each path.
#'             Values closer to 1 will give more cells towards the starting
#'             population, values closer to 0 will give more cells towards the
#'             final population. If a single value is given it will be applied
#'             to all paths.}
#'             \item{\code{[path.nonlinearProb]}}{Probability that a gene
#'             follows a non-linear path along the differentiation path. This
#'             allows more complex gene patterns such as a gene being equally
#'             expressed at the beginning an end of a path but lowly expressed
#'             in the middle.}
#'             \item{\code{[path.sigmaFac]}}{Sigma factor for non-linear gene
#'             paths. A higher value will result in more extreme non-linear
#'             variations along a path.}
#'     }
#'   }
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{splatEstimate}}. For details of the Splat simulation
#' see \code{\link{splatSimulate}}.
#'
#' @name SplatParams
#' @rdname SplatParams
#' @aliases SplatParams-class
#' @exportClass SplatParams
setClass("SplatParams",
         contains = "Params",
         slots = c(nBatches = "numeric",
                   batchCells = "numeric",
                   batch.facLoc = "numeric",
                   batch.facScale = "numeric",
                   batch.rmEffect = "logical",
                   mean.shape = "numeric",
                   mean.rate = "numeric",
                   lib.loc = "numeric",
                   lib.scale = "numeric",
                   lib.norm = "logical",
                   out.prob = "numeric",
                   out.facLoc = "numeric",
                   out.facScale = "numeric",
                   nGroups = "numeric",
                   group.prob = "numeric",
                   de.prob = "numeric",
                   de.downProb = "numeric",
                   de.facLoc = "numeric",
                   de.facScale = "numeric",
                   bcv.common = "numeric",
                   bcv.df = "numeric",
                   dropout.type = "character",
                   dropout.mid = "numeric",
                   dropout.shape = "numeric",
                   path.from = "numeric",
                   path.nSteps = "numeric",
                   path.skew = "numeric",
                   path.nonlinearProb = "numeric",
                   path.sigmaFac = "numeric"),
         prototype = prototype(nBatches = 1,
                               batchCells = 100,
                               batch.facLoc = 0.1,
                               batch.facScale = 0.1,
                               batch.rmEffect = FALSE,
                               mean.rate = 0.3,
                               mean.shape = 0.6,
                               lib.loc = 11,
                               lib.scale = 0.2,
                               lib.norm = FALSE,
                               out.prob = 0.05,
                               out.facLoc = 4,
                               out.facScale = 0.5,
                               nGroups = 1,
                               group.prob = 1,
                               de.prob = 0.1,
                               de.downProb = 0.5,
                               de.facLoc = 0.1,
                               de.facScale = 0.4,
                               bcv.common = 0.1,
                               bcv.df = 60,
                               dropout.type = "none",
                               dropout.mid = 0,
                               dropout.shape = -1,
                               path.from = 0,
                               path.nSteps = 100,
                               path.skew = 0.5,
                               path.nonlinearProb = 0.1,
                               path.sigmaFac = 0.8))

#' The KersplatParams class
#'
#' S4 class that holds parameters for the Kersplat simulation.
#'
#' @section Parameters:
#'
#' The Kersplat simulation uses the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\emph{Mean parameters}}{
#'         \describe{
#'             \item{\code{mean.shape}}{Shape parameter for the mean gamma
#'             distribution.}
#'             \item{\code{mean.rate}}{Rate parameter for the mean gamma
#'             distribution.}
#'             \item{\code{mean.outProb}}{Probability that a gene is an
#'             expression outlier.}
#'             \item{\code{mean.outFacLoc}}{Location (meanlog) parameter for
#'             the expression outlier factor log-normal distribution.}
#'             \item{\code{mean.outFacScale}}{Scale (sdlog) parameter for the
#'             expression outlier factor log-normal distribution.}
#'             \item{\code{mean.dens}}{\code{\link{density}} object describing
#'             the log gene mean density.}
#'             \item{\code{[mean.method]}}{Method to use for simulating gene
#'             means. Either "fit" to sample from a gamma distribution (with
#'             expression outliers) or "density" to sample from the provided
#'             density object.}
#'             \item{\code{[mean.values]}}{Vector of means for each gene.}
#'         }
#'     }
#'     \item{\emph{Biological Coefficient of Variation parameters}}{
#'         \describe{
#'             \item{\code{bcv.common}}{Underlying common dispersion across all
#'             genes.}
#'             \item{\code{[bcv.df]}}{Degrees of Freedom for the BCV inverse
#'             chi-squared distribution.}
#'         }
#'     }
#'     \item{\emph{Network parameters}}{
#'         \describe{
#'             \item{\code{[network.graph]}}{Graph containing the gene network.}
#'             \item{\code{[network.nRegs]}}{Number of regulators in the
#'             network.}
#'         }
#'     }
#'     \item{\emph{Paths parameters}}{
#'         \describe{
#'             \item{\code{[paths.programs]}}{Number of expression programs.}
#'             \item{\code{[paths.design]}}{data.frame describing path
#'             structure. See \code{\link{kersplatSimPaths}} for details.}
#'         }
#'     }
#'     \item{\emph{Library size parameters}}{
#'         \describe{
#'             \item{\code{lib.loc}}{Location (meanlog) parameter for the
#'             library size log-normal distribution, or mean parameter if a
#'             normal distribution is used.}
#'             \item{\code{lib.scale}}{Scale (sdlog) parameter for the library
#'             size log-normal distribution, or sd parameter if a normal
#'             distribution is used.}
#'             \item{\code{lib.dens}}{\code{\link{density}} object describing
#'             the library size density.}
#'             \item{\code{[lib.method]}}{Method to use for simulating library
#'             sizes. Either "fit" to sample from a log-normal distribution or
#'             "density" to sample from the provided density object.}
#'         }
#'     }
#'     \item{\emph{Design parameters}}{
#'         \describe{
#'             \item{\code{[cells.design]}}{data.frame describing cell
#'             structure. See \code{\link{kersplatSimCellMeans}} for details.}
#'         }
#'     }
#'     \item{\emph{Doublet parameters}}{
#'         \describe{
#'             \item{\code{[doublet.prop]}}{Proportion of cells that are
#'             doublets.}
#'         }
#'     }
#'     \item{\emph{Ambient parameters}}{
#'         \describe{
#'             \item{\code{[ambient.scale]}}{Scaling factor for the library
#'             size log-normal distribution when generating ambient library
#'             sizes.}
#'             \item{\code{[ambient.nEmpty]}}{Number of empty cells to
#'             simulate.}
#'         }
#'     }
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{kersplatEstimate}}. For details of the Kersplat simulation
#' see \code{\link{kersplatSimulate}}.
#'
#' @name KersplatParams
#' @rdname KersplatParams
#' @aliases KersplatParams-class
#' @exportClass KersplatParams
setClass("KersplatParams",
         contains = "Params",
         slots = c(mean.shape = "numeric",
                   mean.rate = "numeric",
                   mean.outProb = "numeric",
                   mean.outLoc = "numeric",
                   mean.outScale = "numeric",
                   mean.dens = "density",
                   mean.method = "character",
                   mean.values = "numeric",
                   bcv.common = "numeric",
                   bcv.df = "numeric",
                   network.graph = "ANY",
                   network.nRegs = "numeric",
                   network.regsSet = "logical",
                   paths.nPrograms = "numeric",
                   paths.design = "data.frame",
                   paths.means = "list",
                   lib.loc = "numeric",
                   lib.scale = "numeric",
                   lib.dens = "density",
                   lib.method = "character",
                   cells.design = "data.frame",
                   doublet.prop = "numeric",
                   ambient.scale = "numeric",
                   ambient.nEmpty = "numeric"),
         prototype = prototype(mean.rate = 0.3,
                               mean.shape = 0.6,
                               mean.outProb = 0.05,
                               mean.outLoc = 4,
                               mean.outScale = 0.5,
                               mean.dens = density(rgamma(10000, rate = 0.3,
                                                          shape = 0.6)),
                               mean.method = "fit",
                               mean.values = numeric(),
                               bcv.common = 0.1,
                               bcv.df = 60,
                               network.graph = NULL,
                               network.nRegs = 100,
                               network.regsSet = FALSE,
                               paths.nPrograms = 10,
                               paths.design = data.frame(
                                   Path = 1,
                                   From = 0,
                                   Steps = 100
                               ),
                               paths.means = list(),
                               lib.loc = 11,
                               lib.scale = 0.2,
                               lib.dens = density(rlnorm(10000, 11, 0.2)),
                               lib.method = "fit",
                               cells.design = data.frame(
                                   Path = 1,
                                   Probability = 1,
                                   Alpha = 1,
                                   Beta = 0
                               ),
                               doublet.prop = 0,
                               ambient.scale = 0.05,
                               ambient.nEmpty = 0))

#' The SplatPopParams class
#'
#' S4 class that holds parameters for the splatPop simulation.
#'
#' @section Parameters:
#'
#' In addition to the \code{\link{SplatParams}} parameters, splatPop simulation
#' requires the following parameters:
#'
#' \describe{
#'     \item{\code{[similarity.scale]}}{Scaling factor for pop.cv.param.rate,
#'     where values larger than 1 increase the similarity between individuals in
#'     the population and values less than one make the individuals less
#'     similar.}
#'     \item{\code{[eqtl.n]}}{The number (>1) or percent (<=1) of genes to
#'     assign eQTL effects.}
#'     \item{\code{[eqtl.dist]}}{Maximum distance between eSNP and eGene}
#'     \item{\code{[eqtl.maf.min]}}{Minimum Minor Allele Frequency of eSNPs.}
#'     \item{\code{[eqtl.maf.max]}}{Maximum Minor Allele Frequency of eSNPs.}
#'     \item{\code{[eqtl.group.specific]}}{Percent of eQTL effects to simulate
#'     as group specific.}
#'     \item{\emph{eQTL Effect size distribution parameters. Defaults estimated
#'     from GTEx eQTL mapping results, see vignette for more information.}}{
#'         \describe{
#'             \item{\code{eqtl.ES.shape}}{Shape parameter for the effect size
#'             gamma distribution.}
#'             \item{\code{eqtl.ES.rate}}{Rate parameter for the effect size
#'             gamma distribution.}
#'         }
#'     }
#'     \item{\emph{Bulk Mean Expression distribution parameters. Defaults
#'     estimated from GTEx data, see vignette for more information.}}{
#'         \describe{
#'             \item{\code{pop.mean.shape}}{Shape parameter for the mean (i.e.
#'             bulk) expression gamma distribution}
#'             \item{\code{pop.mean.rate}}{Rate parameter for the mean (i.e.
#'             bulk) expression gamma distribution}
#'         }
#'     }
#'     \item{\emph{Bulk Expression Coefficient of Variation distribution
#'     parameters binned. Defaults estimated from GTEx data, see vignette for
#'     more information.}}{
#'         \describe{
#'             \item{\code{pop.cv.param}}{Dataframe containing gene
#'             mean bin range, and the CV shape, and CV rate parameters for
#'             each of those bins.}
#'         }
#'     }
#'}
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{splatPopEstimate}}. For details of the eQTL simulation
#' see \code{\link{splatPopSimulate}}.
#'
#' @name SplatPopParams
#' @rdname SplatPopParams
#' @aliases SplatPopParams-class
#' @exportClass SplatPopParams
setClass("SplatPopParams",
         contains = "SplatParams",
         slots = c(similarity.scale = "numeric",
                   pop.mean.shape = "numeric",
                   pop.mean.rate = "numeric",
                   pop.cv.bins = "numeric",
                   pop.cv.param = "data.frame",
                   eqtl.n = "numeric",
                   eqtl.dist = "numeric",
                   eqtl.maf.min = "numeric",
                   eqtl.maf.max = "numeric",
                   eqtl.ES.shape = "numeric",
                   eqtl.ES.rate = "numeric",
                   eqtl.group.specific = "numeric"),
         prototype = prototype(similarity.scale = 1.0,
                               pop.mean.shape = 0.3395709,
                               pop.mean.rate = 0.008309486,
                               pop.cv.bins = 10,
                               pop.cv.param =
                                   data.frame(
                                       start = c(0, 0.476, 0.955, 1.86, 3.49,
                                                 6.33, 10.4, 16.3, 26.5,49.9),
                                       end = c(0.476 ,0.955, 1.86, 3.49, 6.33,
                                               10.4, 16.3, 26.5, 49.9, 1e+10),
                                       shape = c(11.636709, 5.084263, 3.161149,
                                                 2.603407, 2.174618, 2.472718,
                                                 2.911565, 3.754947, 3.623545,
                                                 2.540001),
                                       rate = c(8.229737, 3.236401, 1.901426,
                                                1.615142, 1.467896, 2.141105,
                                                3.005807, 4.440894, 4.458207,
                                                2.702462)),
                               eqtl.n = 1,
                               eqtl.dist = 1000000,
                               eqtl.maf.min = 0.05,
                               eqtl.maf.max = 0.5,
                               eqtl.ES.shape = 2.538049,
                               eqtl.ES.rate = 5.962323,
                               eqtl.group.specific = 0.2))

#' The LunParams class
#'
#' S4 class that holds parameters for the Lun simulation.
#'
#' @section Parameters:
#'
#' The Lun simulation uses the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[nGroups]}}{The number of groups to simulate.}
#'     \item{\code{[groupCells]}}{Vector giving the number of cells in each
#'     simulation group/path.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\emph{Mean parameters}}{
#'         \describe{
#'             \item{\code{[mean.shape]}}{Shape parameter for the mean gamma
#'             distribution.}
#'             \item{\code{[mean.rate]}}{Rate parameter for the mean gamma
#'             distribution.}
#'         }
#'     }
#'     \item{\emph{Counts parameters}}{
#'         \describe{
#'             \item{\code{[count.disp]}}{The dispersion parameter for the
#'             counts negative binomial distribution.}
#'         }
#'     }
#'     \item{\emph{Differential expression parameters}}{
#'         \describe{
#'             \item{\code{[de.nGenes]}}{The number of genes that are
#'             differentially expressed in each group}
#'             \item{\code{[de.upProp]}}{The proportion of differentially
#'             expressed genes that are up-regulated in each group}
#'             \item{\code{[de.upFC]}}{The fold change for up-regulated genes}
#'             \item{\code{[de.downFC]}}{The fold change for down-regulated
#'             genes}
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
#' S4 class that holds parameters for the Lun2 simulation.
#'
#' @section Parameters:
#'
#' The Lun2 simulation uses the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\emph{Gene parameters}}{
#'         \describe{
#'             \item{\code{gene.params}}{A \code{data.frame} containing gene
#'             parameters with two columns: \code{Mean} (mean expression for
#'             each gene) and \code{Disp} (dispersion for each gene).}
#'             \item{\code{zi.params}}{A \code{data.frame} containing
#'             zero-inflated gene parameters with three columns: \code{Mean}
#'             (mean expression for each gene), \code{Disp} (dispersion for
#'             each, gene), and \code{Prop} (zero proportion for each gene).}
#'         }
#'     }
#'     \item{\code{[nPlates]}}{The number of plates to simulate.}
#'     \item{\emph{Plate parameters}}{
#'         \describe{
#'             \item{\code{plate.ingroup}}{Character vector giving the plates
#'             considered to be part of the "ingroup".}
#'             \item{\code{plate.mod}}{Plate effect modifier factor. The plate
#'             effect variance is divided by this value.}
#'             \item{\code{plate.var}}{Plate effect variance.}
#'         }
#'     }
#'     \item{\emph{Cell parameters}}{
#'         \describe{
#'             \item{\code{cell.plates}}{Factor giving the plate that each cell
#'             comes from.}
#'             \item{\code{cell.libSizes}}{Library size for each cell.}
#'             \item{\code{cell.libMod}}{Modifier factor for library sizes.
#'             The library sizes are multiplied by this value.}
#'         }
#'     }
#'     \item{\emph{Differential expression parameters}}{
#'         \describe{
#'             \item{\code{de.nGenes}}{Number of differentially expressed
#'             genes.}
#'             \item{\code{de.fc}}{Fold change for differentially expressed
#'             genes.}
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
                   gene.params = "data.frame",
                   zi.params = "data.frame",
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
                               gene.params = data.frame(Mean = rep(3.2, 10000),
                                                        Disp = rep(0.03, 10000)
                                                        ),
                               zi.params = data.frame(Mean = rep(1.6, 10000),
                                                      Disp = rep(0.1, 10000),
                                                      Prop = rep(2.3e-6, 10000)
                                                      ),
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
#'     \item{\code{nGenes}}{The number of genes to simulate (not used).}
#'     \item{\code{nCells}}{The number of cells to simulate in each condition.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\code{SCdat}}{
#'     \code{\link[SingleCellExperiment]{SingleCellExperiment}} containing real
#'     data.}
#'     \item{\code{nDE}}{Number of DE genes to simulate.}
#'     \item{\code{nDP}}{Number of DP genes to simulate.}
#'     \item{\code{nDM}}{Number of DM genes to simulate.}
#'     \item{\code{nDB}}{Number of DB genes to simulate.}
#'     \item{\code{nEE}}{Number of EE genes to simulate.}
#'     \item{\code{nEP}}{Number of EP genes to simulate.}
#'     \item{\code{[sd.range]}}{Interval for fold change standard deviations.}
#'     \item{\code{[modeFC]}}{Values for DP, DM and DB mode fold changes.}
#'     \item{\code{[varInflation]}}{Variance inflation factors for each
#'     condition. If all equal to 1 will be set to \code{NULL} (default).}
#'     \item{\code{[condition]}}{String giving the column that represents
#'     biological group of interest.}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{scDDEstimate}}. See \code{\link[scDD]{simulateSet}} for more
#' details about the parameters. For details of the Splatter implementation of
#' the scDD simulation see \code{\link{scDDSimulate}}.
#'
#' @name SCDDParams
#' @rdname SCDDParams
#' @aliases SCDDParams-class
#' @exportClass SCDDParams
setClass("SCDDParams",
         contains = "Params",
         slots = c(SCdat = "SummarizedExperiment",
                   nDE = "numeric",
                   nDP = "numeric",
                   nDM = "numeric",
                   nDB = "numeric",
                   nEE = "numeric",
                   nEP = "numeric",
                   sd.range = "numeric",
                   modeFC = "numeric",
                   varInflation = "numeric",
                   condition = "character"),
          prototype = prototype(SCdat =
                                   SingleCellExperiment::SingleCellExperiment(),
                               nCells = 100,
                               nDE = 250,
                               nDP = 250,
                               nDM = 250,
                               nDB = 250,
                               nEE = 5000,
                               nEP = 4000,
                               sd.range = c(1, 3),
                               modeFC = c(2, 3, 4),
                               varInflation = c(1, 1),
                               condition = "condition"))

#' The BASiCSParams class
#'
#' S4 class that holds parameters for the BASiCS simulation.
#'
#' @section Parameters:
#'
#' The BASiCS simulation uses the following parameters:
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\emph{Batch parameters}}{
#'         \describe{
#'            \item{\code{nBatches}}{Number of batches to simulate.}
#'            \item{\code{batchCells}}{Number of cells in each batch.}
#'         }
#'     }
#'     \item{\emph{Gene parameters}}{
#'         \describe{
#'             \item{\code{gene.params}}{A \code{data.frame} containing gene
#'             parameters with two columns: \code{Mean} (mean expression for
#'             each biological gene) and \code{Delta} (cell-to-cell
#'             heterogeneity for each biological gene).}
#'         }
#'     }
#'     \item{\emph{Spike-in parameters}}{
#'         \describe{
#'             \item{\code{nSpikes}}{The number of spike-ins to simulate.}
#'             \item{\code{spike.means}}{Input molecules for each spike-in.}
#'         }
#'     }
#'     \item{\emph{Cell parameters}}{
#'         \describe{
#'             \item{\code{cell.params}}{A \code{data.frame} containing gene
#'             parameters with two columns: \code{Phi} (mRNA content factor for
#'             each cell, scaled to sum to the number of cells in each batch)
#'             and \code{S} (capture efficient for each cell).}
#'         }
#'     }
#'     \item{\emph{Variability parameters}}{
#'         \describe{
#'             \item{\code{theta}}{Technical variability parameter for each
#'             batch.}
#'         }
#'     }
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{BASiCSEstimate}}. For details of the BASiCS simulation see
#' \code{\link{BASiCSSimulate}}.
#'
#' @name BASiCSParams
#' @rdname BASiCSParams
#' @aliases BASiCSParams-class
#' @exportClass BASiCSParams
setClass("BASiCSParams",
         contains = "Params",
         slots = c(nBatches = "numeric",
                   batchCells = "numeric",
                   gene.params = "data.frame",
                   nSpikes = "numeric",
                   spike.means = "numeric",
                   cell.params = "data.frame",
                   theta = "numeric"),
         prototype = prototype(nBatches = 1,
                               batchCells = 100,
                               gene.params =
                                   data.frame(
                                       Mean = c(8.36, 10.65, 4.88, 6.29, 21.72,
                                                12.93, 30.19),
                                       Delta = c(1.29, 0.88, 1.51, 1.49, 0.54,
                                                 0.40, 0.85)
                               ),
                               nSpikes = 5,
                               spike.means = c(12.93, 30.19, 1010.72, 7.90,
                                               31.59),
                               cell.params =
                                   data.frame(
                                       Phi = c(1.00, 1.06, 1.09, 1.05, 0.80),
                                       S = c(0.38, 0.40, 0.38, 0.39, 0.34)
                               ),
                               theta = 0.39)
)

#' The MFAParams class
#'
#' S4 class that holds parameters for the mfa simulation.
#'
#' @section Parameters:
#'
#' The mfa simulation uses the following parameters:
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\code{[trans.prop]}}{Proportion of genes that show transient
#'     expression. These genes are briefly up or down-regulated before returning
#'     to their initial state}
#'     \item{\code{[zero.neg]}}{Logical. Whether to set negative expression
#'     values to zero. This will zero-inflate the data.}
#'     \item{\code{[dropout.present]}}{Logical. Whether to simulate dropout.}
#'     \item{\code{dropout.lambda}}{Lambda parameter for the exponential
#'     dropout function.}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{mfaEstimate}}. See \code{\link[mfa]{create_synthetic}} for more
#' details about the parameters. For details of the Splatter implementation of
#' the mfa simulation see \code{\link{mfaSimulate}}.
#'
#' @name MFAParams
#' @rdname MFAParams
#' @aliases MFAParams-class
#' @exportClass MFAParams
setClass("MFAParams",
         contains = "Params",
         slots = c(trans.prop = "numeric",
                   zero.neg = "logical",
                   dropout.present = "logical",
                   dropout.lambda = "numeric"),
         prototype = prototype(trans.prop = 0, zero.neg = TRUE,
                               dropout.present = FALSE, dropout.lambda = 1))


#' The PhenoParams class
#'
#' S4 class that holds parameters for the PhenoPath simulation.
#'
#' @section Parameters:
#'
#' The PhenoPath simulation uses the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\code{[n.de]}}{Number of genes to simulate from the differential
#'     expression regime}
#'     \item{\code{[n.pst]}}{Number of genes to simulate from the pseudotime
#'     regime}
#'     \item{\code{[n.pst.beta]}}{Number of genes to simulate from the
#'     pseudotime + beta interactions regime}
#'     \item{\code{[n.de.pst.beta]}}{Number of genes to simulate from the
#'     differential expression + pseudotime + interactions regime}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{phenoEstimate}}. For details of the PhenoPath simulation
#' see \code{\link{phenoSimulate}}.
#'
#' @name PhenoParams
#' @rdname PhenoParams
#' @aliases PhenoParams-class
#' @exportClass PhenoParams
setClass("PhenoParams",
         contains = "Params",
         slots = c(n.de = "numeric",
                   n.pst = "numeric",
                   n.pst.beta = "numeric",
                   n.de.pst.beta = "numeric"),
         prototype = prototype(n.de = 2500, n.pst = 2500, n.pst.beta = 2500,
                               n.de.pst.beta = 2500))


#' The ZINBParams class
#'
#' S4 class that holds parameters for the ZINB-WaVE simulation.
#'
#' @section Parameters:
#'
#' The ZINB-WaVE simulation uses the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\code{model}}{Object describing a ZINB model.}
#' }
#'
#' The majority of the parameters for this simulation are stored in a
#' \code{\link[zinbwave]{ZinbModel}} object. Please refer to the documentation
#' for this class and its constructor(\code{\link[zinbwave]{zinbModel}}) for
#' details about all the parameters.
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{zinbEstimate}}. For details of the ZINB-WaVE simulation
#' see \code{\link{zinbSimulate}}.
#'
#' @name ZINBParams
#' @rdname ZINBParams
#' @aliases ZINBParams-class
#' @exportClass ZINBParams
setClass("ZINBParams",
         contains = "Params",
         slots = c(model = "ANY"),
         prototype = prototype(nGenes = 100, nCells = 50))


#' The SparseDCParams class
#'
#' S4 class that holds parameters for the SparseDC simulation.
#'
#' @section Parameters:
#'
#' The SparseDC simulation uses the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate in each condition.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\code{markers.n}}{Number of marker genes to simulate for each
#'     cluster.}
#'     \item{\code{markers.shared}}{Number of marker genes for each cluster
#'     shared between conditions. Must be less than or equal to
#'     \code{markers.n}}.
#'     \item{\code{[markers.same]}}{Logical. Whether each cluster should have
#'     the same set of marker genes.}
#'     \item{\code{clusts.c1}}{Numeric vector of clusters present in
#'     condition 1. The number of times a cluster is repeated controls the
#'     proportion of cells from that cluster.}
#'     \item{\code{clusts.c2}}{Numeric vector of clusters present in
#'     condition 2. The number of times a cluster is repeated controls the
#'     proportion of cells from that cluster.}
#'     \item{\code{[mean.lower]}}{Lower bound for cluster gene means.}
#'     \item{\code{[mean.upper]}}{Upper bound for cluster gene means.}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{sparseDCEstimate}}. For details of the SparseDC simulation
#' see \code{\link{sparseDCSimulate}}.
#'
#' @name SparseDCParams
#' @rdname SparseDCParams
#' @aliases SparseDCParams-class
#' @exportClass SparseDCParams
setClass("SparseDCParams",
         contains = "Params",
         slots = c(markers.n = "numeric",
                   markers.shared = "numeric",
                   markers.same = "logical",
                   clusts.c1 = "numeric",
                   clusts.c2 = "numeric",
                   mean.lower = "numeric",
                   mean.upper = "numeric"),
         prototype = prototype(markers.n = 0,
                               markers.shared = 0,
                               markers.same = FALSE,
                               clusts.c1 = 1,
                               clusts.c2 = 1,
                               mean.lower = 1,
                               mean.upper = 2))
