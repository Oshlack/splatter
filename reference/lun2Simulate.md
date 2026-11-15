# Lun2 simulation

Simulate single-cell RNA-seq count data using the method described in
Lun and Marioni "Overcoming confounding plate effects in differential
expression analyses of single-cell RNA-seq data".

## Usage

``` r
lun2Simulate(
  params = newLun2Params(),
  zinb = FALSE,
  sparsify = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- params:

  Lun2Params object containing simulation parameters.

- zinb:

  logical. Whether to use a zero-inflated model.

- sparsify:

  logical. Whether to automatically convert assays to sparse matrices if
  there will be a size reduction.

- verbose:

  logical. Whether to print progress messages.

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

SingleCellExperiment containing simulated counts.

## Details

The Lun2 simulation uses a negative-binomial distribution where the
means and dispersions have been sampled from a real dataset (using
[`lun2Estimate`](http://oshlacklab.com/splatter/reference/lun2Estimate.md)).
The other core feature of the Lun2 simulation is the addition of plate
effects. Differential expression can be added between two groups of
plates (an "ingroup" and all other plates). Library size factors are
also applied and optionally a zero-inflated negative-binomial can be
used.

If the number of genes to simulate differs from the number of provided
gene parameters or the number of cells to simulate differs from the
number of library sizes the relevant parameters will be sampled with a
warning. This allows any number of genes or cells to be simulated
regardless of the number in the dataset used in the estimation step but
has the downside that some genes or cells may be simulated multiple
times.

## References

Lun ATL, Marioni JC. Overcoming confounding plate effects in
differential expression analyses of single-cell RNA-seq data.
Biostatistics (2017).

Paper:
[dx.doi.org/10.1093/biostatistics/kxw055](http://oshlacklab.com/splatter/reference/dx.doi.org/10.1093/biostatistics/kxw055)

Code: <https://github.com/MarioniLab/PlateEffects2016>

## Examples

``` r
sim <- lun2Simulate()
#> Getting parameters...
#> Simulating plate means...
#> Simulating library size factors...
#> Simulating cell means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.13 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.5 * dense matrix
#> Converting 'TrueCounts' to sparse matrix: estimated sparse size 0.13 * dense matrix
#> Done!
```
