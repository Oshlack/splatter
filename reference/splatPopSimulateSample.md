# splatPopSimulateSample simulation

Simulate count data for one sample from a fictional single-cell RNA-seq
experiment using the Splat method.

## Usage

``` r
splatPopSimulateSample(
  params = newSplatPopParams(),
  method = c("single", "groups", "paths"),
  batch = "batch1",
  counts.only = FALSE,
  verbose = TRUE,
  sample.means,
  ...
)
```

## Arguments

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- method:

  which simulation method to use. Options are "single" which produces a
  single population, "groups" which produces distinct groups (eg. cell
  types), "paths" which selects cells from continuous trajectories (eg.
  differentiation processes).

- batch:

  Batch number.

- counts.only:

  logical. Whether to return only the counts.

- verbose:

  logical. Whether to print progress messages.

- sample.means:

  Gene means to use if running splatSimulatePop().

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

SingleCellExperiment object containing the simulated counts and
intermediate values for one sample.

## Details

This function closely mirrors
[`splatSimulate`](http://oshlacklab.com/splatter/reference/splatSimulate.md).
The main difference is that it takes the means simulated by
splatPopSimulateMeans instead of randomly sampling a mean for each gene.
For details about this function see the documentation for
[`splatSimulate`](http://oshlacklab.com/splatter/reference/splatSimulate.md).

## See also

[`splatSimLibSizes`](http://oshlacklab.com/splatter/reference/splatSimLibSizes.md),
[`splatPopSimGeneMeans`](http://oshlacklab.com/splatter/reference/splatPopSimGeneMeans.md),
[`splatSimBatchEffects`](http://oshlacklab.com/splatter/reference/splatSimBatchEffects.md),
[`splatSimBatchCellMeans`](http://oshlacklab.com/splatter/reference/splatSimBatchCellMeans.md),
[`splatSimDE`](http://oshlacklab.com/splatter/reference/splatSimDE.md),
[`splatSimCellMeans`](http://oshlacklab.com/splatter/reference/splatSimCellMeans.md),
[`splatSimBCVMeans`](http://oshlacklab.com/splatter/reference/splatSimBCVMeans.md),
[`splatSimTrueCounts`](http://oshlacklab.com/splatter/reference/splatSimTrueCounts.md),
[`splatSimDropout`](http://oshlacklab.com/splatter/reference/splatSimDropout.md),
[`splatPopSimulateSC`](http://oshlacklab.com/splatter/reference/splatPopSimulateSC.md)
