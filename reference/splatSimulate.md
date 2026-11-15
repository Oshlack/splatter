# Splat simulation

Simulate count data from a fictional single-cell RNA-seq experiment
using the Splat method.

## Usage

``` r
splatSimulate(
  params = newSplatParams(),
  method = c("single", "groups", "paths"),
  sparsify = TRUE,
  verbose = TRUE,
  ...
)

splatSimulateSingle(params = newSplatParams(), verbose = TRUE, ...)

splatSimulateGroups(params = newSplatParams(), verbose = TRUE, ...)

splatSimulatePaths(params = newSplatParams(), verbose = TRUE, ...)
```

## Arguments

- params:

  SplatParams object containing parameters for the simulation. See
  [`SplatParams`](http://oshlacklab.com/splatter/reference/SplatParams.md)
  for details.

- method:

  which simulation method to use. Options are "single" which produces a
  single population, "groups" which produces distinct groups (eg. cell
  types), or "paths" which selects cells from continuous trajectories
  (eg. differentiation processes).

- sparsify:

  logical. Whether to automatically convert assays to sparse matrices if
  there will be a size reduction.

- verbose:

  logical. Whether to print progress messages.

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

SingleCellExperiment object containing the simulated counts and
intermediate values.

## Details

Parameters can be set in a variety of ways. If no parameters are
provided the default parameters are used. Any parameters in `params` can
be overridden by supplying additional arguments through a call to
[`setParams`](http://oshlacklab.com/splatter/reference/setParams.md).
This design allows the user flexibility in how they supply parameters
and allows small adjustments without creating a new `SplatParams`
object. See examples for a demonstration of how this can be used.

The simulation involves the following steps:

1.  Set up simulation object

2.  Simulate library sizes

3.  Simulate gene means

4.  Simulate groups/paths

5.  Simulate BCV adjusted cell means

6.  Simulate true counts

7.  Simulate dropout

8.  Create final dataset

The final output is a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object that contains the simulated counts but also the values for
various intermediate steps. These are stored in the
[`colData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
(for cell specific information),
[`rowData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
(for gene specific information) or
[`assays`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
(for gene by cell matrices) slots. This additional information includes:

- `colData`:

  Cell

  : Unique cell identifier.

  Group

  : The group or path the cell belongs to.

  ExpLibSize

  : The expected library size for that cell.

  Step (paths only)

  : how far along the path each cell is.

- `rowData`:

  Gene

  : Unique gene identifier.

  BaseGeneMean

  : The base expression level for that gene.

  OutlierFactor

  : Expression outlier factor for that gene. Values of 1 indicate the
    gene is not an expression outlier.

  GeneMean

  : Expression level after applying outlier factors.

  BatchFac\[Batch\]

  : The batch effects factor for each gene for a particular batch.

  DEFac\[Group\]

  : The differential expression factor for each gene in a particular
    group. Values of 1 indicate the gene is not differentially
    expressed.

  SigmaFac\[Path\]

  : Factor applied to genes that have non-linear changes in expression
    along a path.

- `assays`:

  BatchCellMeans

  : The mean expression of genes in each cell after adding batch
    effects.

  BaseCellMeans

  : The mean expression of genes in each cell after any differential
    expression and adjusted for expected library size.

  BCV

  : The Biological Coefficient of Variation for each gene in each cell.

  CellMeans

  : The mean expression level of genes in each cell adjusted for BCV.

  TrueCounts

  : The simulated counts before dropout.

  Dropout

  : Logical matrix showing which values have been dropped in which
    cells.

Values that have been added by Splatter are named using `UpperCamelCase`
in order to differentiate them from the values added by analysis
packages which typically use `underscore_naming`.

## References

Zappia L, Phipson B, Oshlack A. Splatter: simulation of single-cell RNA
sequencing data. Genome Biology (2017).

Paper:
[10.1186/s13059-017-1305-0](http://oshlacklab.com/splatter/reference/10.1186/s13059-017-1305-0)

Code: <https://github.com/Oshlack/splatter>

## See also

[`splatSimLibSizes`](http://oshlacklab.com/splatter/reference/splatSimLibSizes.md),
[`splatSimGeneMeans`](http://oshlacklab.com/splatter/reference/splatSimGeneMeans.md),
[`splatSimBatchEffects`](http://oshlacklab.com/splatter/reference/splatSimBatchEffects.md),
[`splatSimBatchCellMeans`](http://oshlacklab.com/splatter/reference/splatSimBatchCellMeans.md),
[`splatSimDE`](http://oshlacklab.com/splatter/reference/splatSimDE.md),
[`splatSimCellMeans`](http://oshlacklab.com/splatter/reference/splatSimCellMeans.md),
[`splatSimBCVMeans`](http://oshlacklab.com/splatter/reference/splatSimBCVMeans.md),
[`splatSimTrueCounts`](http://oshlacklab.com/splatter/reference/splatSimTrueCounts.md),
[`splatSimDropout`](http://oshlacklab.com/splatter/reference/splatSimDropout.md)

## Examples

``` r
# Simulation with default parameters
sim <- splatSimulate()
#> Getting parameters...
#> Creating simulation object...
#> Simulating library sizes...
#> Simulating gene means...
#> Simulating BCV...
#> Simulating counts...
#> Simulating dropout (if needed)...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'BatchCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BaseCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BCV': estimated sparse size 1.5 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.49 * dense matrix
#> Skipping 'TrueCounts': estimated sparse size 1.63 * dense matrix
#> Skipping 'counts': estimated sparse size 1.63 * dense matrix
#> Done!

# \donttest{
# Simulation with different number of genes
sim <- splatSimulate(nGenes = 1000)
#> Getting parameters...
#> Creating simulation object...
#> Simulating library sizes...
#> Simulating gene means...
#> Simulating BCV...
#> Simulating counts...
#> Simulating dropout (if needed)...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'BatchCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BaseCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BCV': estimated sparse size 1.5 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'TrueCounts': estimated sparse size 2.58 * dense matrix
#> Skipping 'counts': estimated sparse size 2.58 * dense matrix
#> Done!
# Simulation with custom parameters
params <- newSplatParams(nGenes = 100, mean.rate = 0.5)
sim <- splatSimulate(params)
#> Getting parameters...
#> Creating simulation object...
#> Simulating library sizes...
#> Simulating gene means...
#> Simulating BCV...
#> Simulating counts...
#> Simulating dropout (if needed)...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'BatchCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BaseCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BCV': estimated sparse size 1.5 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'TrueCounts': estimated sparse size 2.93 * dense matrix
#> Skipping 'counts': estimated sparse size 2.93 * dense matrix
#> Done!
# Simulation with adjusted custom parameters
sim <- splatSimulate(params, mean.rate = 0.6, out.prob = 0.2)
#> Getting parameters...
#> Creating simulation object...
#> Simulating library sizes...
#> Simulating gene means...
#> Simulating BCV...
#> Simulating counts...
#> Simulating dropout (if needed)...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'BatchCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BaseCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BCV': estimated sparse size 1.5 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'TrueCounts': estimated sparse size 2.85 * dense matrix
#> Skipping 'counts': estimated sparse size 2.85 * dense matrix
#> Done!
# Simulate groups
sim <- splatSimulate(method = "groups")
#> Getting parameters...
#> Warning: nGroups is 1, switching to single mode
#> Creating simulation object...
#> Simulating library sizes...
#> Simulating gene means...
#> Simulating BCV...
#> Simulating counts...
#> Simulating dropout (if needed)...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'BatchCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BaseCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BCV': estimated sparse size 1.5 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.49 * dense matrix
#> Skipping 'TrueCounts': estimated sparse size 1.66 * dense matrix
#> Skipping 'counts': estimated sparse size 1.66 * dense matrix
#> Done!
# Simulate paths
sim <- splatSimulate(method = "paths")
#> Getting parameters...
#> Creating simulation object...
#> Simulating library sizes...
#> Simulating gene means...
#> Simulating path endpoints...
#> Simulating path steps...
#> Simulating BCV...
#> Simulating counts...
#> Simulating dropout (if needed)...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'BatchCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BaseCellMeans': estimated sparse size 1.5 * dense matrix
#> Skipping 'BCV': estimated sparse size 1.5 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.49 * dense matrix
#> Skipping 'TrueCounts': estimated sparse size 1.65 * dense matrix
#> Skipping 'counts': estimated sparse size 1.65 * dense matrix
#> Done!
# }
```
