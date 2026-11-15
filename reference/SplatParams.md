# The SplatParams class

S4 class that holds parameters for the Splat simulation.

## Parameters

The Splat simulation requires the following parameters:

- `nGenes`:

  The number of genes to simulate.

- `nCells`:

  The number of cells to simulate.

- `[seed]`:

  Seed to use for generating random numbers.

- *Batch parameters*:

  `[nBatches]`

  : The number of batches to simulate.

  `[batchCells]`

  : Vector giving the number of cells in each batch.

  `[batch.facLoc]`

  : Location (meanlog) parameter for the batch effect factor log-normal
    distribution. Can be a vector.

  `[batch.facScale]`

  : Scale (sdlog) parameter for the batch effect factor log-normal
    distribution. Can be a vector.

  `[batch.rmEffect]`

  : Logical, removes the batch effect and continues with the simulation
    when TRUE. This allows the user to test batch removal algorithms
    without having to calculate the new expected cell means with batch
    removed.

- *Mean parameters*:

  `mean.shape`

  : Shape parameter for the mean gamma distribution.

  `mean.rate`

  : Rate parameter for the mean gamma distribution.

- *Library size parameters*:

  `lib.loc`

  : Location (meanlog) parameter for the library size log-normal
    distribution, or mean parameter if a normal distribution is used.

  `lib.scale`

  : Scale (sdlog) parameter for the library size log-normal
    distribution, or sd parameter if a normal distribution is used.

  `lib.norm`

  : Logical. Whether to use a normal distribution for library sizes
    instead of a log-normal.

- *Expression outlier parameters*:

  `out.prob`

  : Probability that a gene is an expression outlier.

  `out.facLoc`

  : Location (meanlog) parameter for the expression outlier factor
    log-normal distribution.

  `out.facScale`

  : Scale (sdlog) parameter for the expression outlier factor log-normal
    distribution.

- *Group parameters*:

  `[nGroups]`

  : The number of groups or paths to simulate.

  `[group.prob]`

  : Probability that a cell comes from a group.

- *Differential expression parameters*:

  `[de.prob]`

  : Probability that a gene is differentially expressed in a group. Can
    be a vector.

  `[de.downProb]`

  : Probability that a differentially expressed gene is down-regulated.
    Can be a vector.

  `[de.facLoc]`

  : Location (meanlog) parameter for the differential expression factor
    log-normal distribution. Can be a vector.

  `[de.facScale]`

  : Scale (sdlog) parameter for the differential expression factor
    log-normal distribution. Can be a vector.

- *Biological Coefficient of Variation parameters*:

  `bcv.common`

  : Underlying common dispersion across all genes.

  `bcv.df`

  : Degrees of Freedom for the BCV inverse chi-squared distribution.

- *Dropout parameters*:

  `dropout.type`

  : The type of dropout to simulate. "none" indicates no dropout,
    "experiment" is global dropout using the same parameters for every
    cell, "batch" uses the same parameters for every cell in each batch,
    "group" uses the same parameters for every cell in each groups and
    "cell" uses a different set of parameters for each cell.

  `dropout.mid`

  : Midpoint parameter for the dropout logistic function.

  `dropout.shape`

  : Shape parameter for the dropout logistic function.

- *Differentiation path parameters*:

  `[path.from]`

  : Vector giving the originating point of each path. This allows path
    structure such as a cell type which differentiates into an
    intermediate cell type that then differentiates into two mature cell
    types. A path structure of this form would have a "from" parameter
    of c(0, 1, 1) (where 0 is the origin). If no vector is given all
    paths will start at the origin.

  `[path.nSteps]`

  : Vector giving the number of steps to simulate along each path. If a
    single value is given it will be applied to all paths. This
    parameter was previously called `path.length`.

  `[path.skew]`

  : Vector giving the skew of each path. Values closer to 1 will give
    more cells towards the starting population, values closer to 0 will
    give more cells towards the final population. If a single value is
    given it will be applied to all paths.

  `[path.nonlinearProb]`

  : Probability that a gene follows a non-linear path along the
    differentiation path. This allows more complex gene patterns such as
    a gene being equally expressed at the beginning an end of a path but
    lowly expressed in the middle.

  `[path.sigmaFac]`

  : Sigma factor for non-linear gene paths. A higher value will result
    in more extreme non-linear variations along a path.

The parameters not shown in brackets can be estimated from real data
using
[`splatEstimate`](http://oshlacklab.com/splatter/reference/splatEstimate.md).
For details of the Splat simulation see
[`splatSimulate`](http://oshlacklab.com/splatter/reference/splatSimulate.md).
