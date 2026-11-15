# Kersplat sample

Sample cells for the Kersplat simulation

## Usage

``` r
kersplatSample(params, sparsify = TRUE, verbose = TRUE)
```

## Arguments

- params:

  KersplatParams object containing simulation parameters.

- sparsify:

  logical. Whether to automatically convert assays to sparse matrices if
  there will be a size reduction.

- verbose:

  logical. Whether to print progress messages

## Value

SingleCellExperiment object containing the simulated counts and
intermediate values.

## Details

The second stage is a two-step Kersplat simulation is to generate cells
based on a complete
[`KersplatParams`](http://oshlacklab.com/splatter/reference/KersplatParams.md)
object. intermediate parameters.

The sampling process involves the following steps:

1.  Simulate library sizes for each cell

2.  Simulate means for each cell

3.  Simulate endogenous counts for each cell

4.  Simulate ambient counts for each cell

5.  Simulate final counts for each cell

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

  Type

  : Whether the cell is a Cell, Doublet or Empty.

  CellLibSize

  : The expected number of endogenous counts for that cell.

  AmbientLibSize

  : The expected number of ambient counts for that cell.

  Path

  : The path the cell belongs to.

  Step

  : How far along the path each cell is.

  Path1

  : For doublets the path of the first partner in the doublet (otherwise
    `NA`).

  Step1

  : For doublets the step of the first partner in the doublet (otherwise
    `NA`).

  Path2

  : For doublets the path of the second partner in the doublet
    (otherwise `NA`).

  Step2

  : For doublets the step of the second partner in the doublet
    (otherwise `NA`).

- `rowData`:

  Gene

  : Unique gene identifier.

  BaseMean

  : The base expression level for that gene.

  AmbientMean

  : The ambient expression level for that gene.

- `assays`:

  CellMeans

  : The mean expression of genes in each cell after any differential
    expression and adjusted for expected library size.

  CellCounts

  : Endogenous count matrix.

  AmbientCounts

  : Ambient count matrix.

  counts

  : Final count matrix.

Values that have been added by Splatter are named using `UpperCamelCase`
in order to differentiate them from the values added by analysis
packages which typically use `underscore_naming`.

## See also

[`kersplatSimLibSizes`](http://oshlacklab.com/splatter/reference/kersplatSimLibSizes.md),
[`kersplatSimCellMeans`](http://oshlacklab.com/splatter/reference/kersplatSimCellMeans.md),
[`kersplatSimCellCounts`](http://oshlacklab.com/splatter/reference/kersplatSimCellCounts.md),
[`kersplatSimAmbientCounts`](http://oshlacklab.com/splatter/reference/kersplatSimAmbientCounts.md),
[`kersplatSimCounts`](http://oshlacklab.com/splatter/reference/kersplatSimCounts.md)

## Examples

``` r

if (requireNamespace("igraph", quietly = TRUE)) {
    params <- kersplatSetup()
    sim <- kersplatSample(params)
}
#> Setting up parameters...
#> Generating gene network...
#> Selecting regulators...
#> Simulating means...
#> Sampling from gamma distribution...
#> Simulating paths...
#> Simulating path 1...
#> Creating simulation object...
#> Simulating library sizes...
#> Sampling from log-normal distribution...
#> Assigning cells to paths...
#> Assigning cells to steps...
#> Simulating cell means...
#> Applying BCV adjustment...
#> Simulating cell counts...
#> Simulating ambient counts...
#> Simulating final counts...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'CellMeans': estimated sparse size 1.49 * dense matrix
#> Skipping 'CellCounts': estimated sparse size 1.65 * dense matrix
#> Converting 'AmbientCounts' to sparse matrix: estimated sparse size 0.43 * dense matrix
#> Skipping 'counts': estimated sparse size 1.65 * dense matrix
```
