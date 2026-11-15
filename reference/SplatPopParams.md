# The SplatPopParams class

S4 class that holds parameters for the splatPop simulation.

## Parameters

In addition to the
[`SplatParams`](http://oshlacklab.com/splatter/reference/SplatParams.md)
parameters, splatPop simulation requires the following parameters:

- `[similarity.scale]`:

  Scaling factor for pop.cv.param.rate, where values larger than 1
  increase the similarity between individuals in the population and
  values less than one make the individuals less similar.

- `[eqtl.n]`:

  The number (\>1) or percent (\<=1) of genes to assign eQTL effects.

- `[eqtl.dist]`:

  Maximum distance between eSNP and eGene

- `[eqtl.maf.min]`:

  Minimum Minor Allele Frequency of eSNPs.

- `[eqtl.maf.max]`:

  Maximum Minor Allele Frequency of eSNPs.

- `[eqtl.coreg]`:

  Proportion of eGenes to have a shared eSNP (i.e., co-regulated genes)

- `[eqtl.group.specific]`:

  Percent of eQTL effects to simulate as group specific.

- `[eqtl.condition.specific]`:

  Percent of eQTL effects to simulate as condition specific.

- *eQTL Effect size distribution parameters. Defaults estimated from
  GTEx eQTL mapping results, see vignette for more information.*:

  `eqtl.ES.shape`

  : Shape parameter for the effect size gamma distribution.

  `eqtl.ES.rate`

  : Rate parameter for the effect size gamma distribution.

- *Bulk Mean Expression distribution parameters. Defaults estimated from
  GTEx data, see vignette for more information.*:

  `pop.mean.shape`

  : Shape parameter for the mean (i.e. bulk) expression gamma
    distribution

  `pop.mean.rate`

  : Rate parameter for the mean (i.e. bulk) expression gamma
    distribution

- *Bulk Expression Coefficient of Variation distribution parameters
  binned. Defaults estimated from GTEx data, see vignette for more
  information.*:

  `pop.cv.param`

  : Dataframe containing gene mean bin range, and the CV shape, and CV
    rate parameters for each of those bins.

- *Specify number of samples per batch. Note that splatPop will randomly
  assign donors to be present in multiple batches to fulfill the
  specified nBatches and batch.size parameters. For example, if 10
  samples are simulated with batchPool.n=4 and batchPool.size= 4, then 6
  samples will be randomly chosen to be replicated in two pools.*:

  `batch.size`

  : The number of donors in each pool/batch.

- *Specify shape and rate of gamma distribution to sample number of
  cells per batch per donor. Will only be used if nCells parameter is
  set to 0.*:

  `nCells.sample`

  : True/False if nCells should be set as nCells or sampled from a gamma
    distribution for each batch/donor.

  `nCells.shape`

  : Shape parameter for the nCells per batch per donor distribution.

  `nCells.rate`

  : Rate parameter for the nCells per batch per donor distribution.

- *Condition/treatment differential expression parameters*:

  `[nConditions]`

  : The number of conditions/treatments to divide samples into.

  `[condition.prob]`

  : Probability that a sample belongs to each condition/treatment group.
    Can be a vector.

  `[cde.prob]`

  : Probability that a gene is differentially expressed in a condition
    group. Can be a vector.

  `[cde.downProb]`

  : Probability that a conditionally differentially expressed gene is
    down-regulated. Can be a vector.

  `[cde.facLoc]`

  : Location (meanlog) parameter for the conditional differential
    expression factor log-normal distribution. Can be a vector.

  `[cde.facScale]`

  : Scale (sdlog) parameter for the conditional differential expression
    factor log-normal distribution. Can be a vector.

The parameters not shown in brackets can be estimated from real data
using
[`splatPopEstimate`](http://oshlacklab.com/splatter/reference/splatPopEstimate.md).
For details of the eQTL simulation see
[`splatPopSimulate`](http://oshlacklab.com/splatter/reference/splatPopSimulate.md).
