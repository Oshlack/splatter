# BASiCS simulation

Simulate counts using the BASiCS method.

## Usage

``` r
BASiCSSimulate(
  params = newBASiCSParams(),
  sparsify = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- params:

  BASiCSParams object containing simulation parameters.

- sparsify:

  logical. Whether to automatically convert assays to sparse matrices if
  there will be a size reduction.

- verbose:

  logical. Whether to print progress messages

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

SingleCellExperiment containing simulated counts

## Details

This function is just a wrapper around
[`BASiCS_Sim`](https://rdrr.io/pkg/BASiCS/man/BASiCS_Sim.html) that
takes a
[`BASiCSParams`](http://oshlacklab.com/splatter/reference/BASiCSParams.md),
runs the simulation then converts the output to a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object. See
[`BASiCS_Sim`](https://rdrr.io/pkg/BASiCS/man/BASiCS_Sim.html) for more
details of how the simulation works.

## References

Vallejos CA, Marioni JC, Richardson S. BASiCS: Bayesian Analysis of
Single-Cell Sequencing data. PLoS Computational Biology (2015).

Paper:
[10.1371/journal.pcbi.1004333](http://oshlacklab.com/splatter/reference/10.1371/journal.pcbi.1004333)

Code: <https://github.com/catavallejos/BASiCS>

## Examples

``` r
if (requireNamespace("BASiCS", quietly = TRUE)) {
    sim <- BASiCSSimulate()
}
#> Getting parameters...
#> Warning: Number of gene.params not equal to nGenes, gene.params will be sampled.
#> Warning: Number of cell.params not equal to nCells, cell.params will be sampled.
#> Simulating counts with BASiCS...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Skipping 'counts': estimated sparse size 1.08 * dense matrix
#> Done!
```
