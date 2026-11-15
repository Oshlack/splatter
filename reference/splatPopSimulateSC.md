# splatPopSimulateSC

Simulate count data for a population from a fictional single-cell
RNA-seq experiment using the Splat method.

## Usage

``` r
splatPopSimulateSC(
  sim.means,
  params,
  key,
  method = c("single", "groups", "paths"),
  counts.only = FALSE,
  conditions = NULL,
  sparsify = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- sim.means:

  Matrix or list of matrices of gene means for the population. Output
  from
  [`splatPopSimulateMeans()`](http://oshlacklab.com/splatter/reference/splatPopSimulateMeans.md).

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- key:

  data.frame object containing a full or partial splatPop key. Output
  from
  [`splatPopSimulateMeans()`](http://oshlacklab.com/splatter/reference/splatPopSimulateMeans.md).

- method:

  which simulation method to use. Options are "single" which produces a
  single cell population for each sample, "groups" which produces
  distinct groups (eg. cell types) for each sample (note, this creates
  separate groups from those created in `popSimulate` with only DE
  effects), and "paths" which selects cells from continuous trajectories
  (eg. differentiation processes).

- counts.only:

  logical. Whether to return only the counts.

- conditions:

  named array with conditional group assignment for each sample. Output
  from
  [`splatPopSimulateMeans()`](http://oshlacklab.com/splatter/reference/splatPopSimulateMeans.md).

- sparsify:

  logical. Whether to automatically convert assays to sparse matrices if
  there will be a size reduction.

- verbose:

  logical. Whether to print progress messages.

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

SingleCellExperiment object containing simulated counts, intermediate
values like the gene means simulated in `splatPopSimulateMeans`, and
information about the differential expression and eQTL effects assigned
to each gene.

## Examples

``` r
# \donttest{
if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
    requireNamespace("preprocessCore", quietly = TRUE)) {
    params <- newSplatPopParams()
    sim.means <- splatPopSimulateMeans()
    sim <- splatPopSimulateSC(sim.means$means, params, sim.means$key)
}
#> Simulating data for genes in GFF...
#> Simulating gene means for population...
#> Simulating sc counts for Group1...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'BatchCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BaseCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BCV': estimated sparse size 1.5 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.49 * dense matrix
#> Skipping 'TrueCounts': estimated sparse size 2.74 * dense matrix
#> Skipping 'counts': estimated sparse size 2.74 * dense matrix
#> Done!
# }
```
