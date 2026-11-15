# scDD simulation

Simulate counts using the scDD method.

## Usage

``` r
scDDSimulate(
  params = newSCDDParams(),
  plots = FALSE,
  plot.file = NULL,
  sparsify = TRUE,
  verbose = TRUE,
  BPPARAM = SerialParam(),
  ...
)
```

## Arguments

- params:

  SCDDParams object containing simulation parameters.

- plots:

  logical. whether to generate scDD fold change and validation plots.

- plot.file:

  File path to save plots as PDF.

- sparsify:

  logical. Whether to automatically convert assays to sparse matrices if
  there will be a size reduction.

- verbose:

  logical. Whether to print progress messages

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance giving the parallel back-end to be used. Default is
  [`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  which uses a single core.

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

SingleCellExperiment containing simulated counts

## Details

This function is just a wrapper around
[`simulateSet`](https://rdrr.io/pkg/scDD/man/simulateSet.html) that
takes a
[`SCDDParams`](http://oshlacklab.com/splatter/reference/SCDDParams.md),
runs the simulation then converts the output to a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object. See
[`simulateSet`](https://rdrr.io/pkg/scDD/man/simulateSet.html) for more
details about how the simulation works.

## References

Korthauer KD, Chu L-F, Newton MA, Li Y, Thomson J, Stewart R, et al. A
statistical approach for identifying differential distributions in
single-cell RNA-seq experiments. Genome Biology (2016).

Paper:
[10.1186/s13059-016-1077-y](http://oshlacklab.com/splatter/reference/10.1186/s13059-016-1077-y)

Code: <https://github.com/kdkorthauer/scDD>

## Examples

``` r
# \donttest{
sim <- scDDSimulate()
#> Simulating counts with scDD...
#> Setting up parallel back-end using 1 cores
#> Identifying a set of genes to simulate from...
#> Simulating DE fold changes...
#> Simulating individual genes...
#> Done! Simulated 250 DE, 250 DP, 250 DM, 250 DB, 5000 EE, and 4000 EP genes 
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'counts': estimated sparse size 1.36 * dense matrix
#> Done!
# }
```
