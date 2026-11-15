# Lun simulation

Simulate single-cell RNA-seq count data using the method described in
Lun, Bach and Marioni "Pooling across cells to normalize single-cell RNA
sequencing data with many zero counts".

## Usage

``` r
lunSimulate(params = newLunParams(), sparsify = TRUE, verbose = TRUE, ...)
```

## Arguments

- params:

  LunParams object containing Lun simulation parameters.

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

The Lun simulation generates gene mean expression levels from a gamma
distribution with `shape = mean.shape` and `rate = mean.rate`. Counts
are then simulated from a negative binomial distribution with
`mu = means` and `size = 1 / bcv.common`. In addition each cell is given
a size factor (`2 ^ rnorm(nCells, mean = 0, sd = 0.5)`) and differential
expression can be simulated with fixed fold changes.

See [`LunParams`](http://oshlacklab.com/splatter/reference/LunParams.md)
for details of the parameters.

## References

Lun ATL, Bach K, Marioni JC. Pooling across cells to normalize
single-cell RNA sequencing data with many zero counts. Genome Biology
(2016).

Paper:
[dx.doi.org/10.1186/s13059-016-0947-7](http://oshlacklab.com/splatter/reference/dx.doi.org/10.1186/s13059-016-0947-7)

Code: <https://github.com/MarioniLab/Deconvolution2016>

## Examples

``` r
sim <- lunSimulate()
#> Getting parameters...
#> Simulating means...
#> Simulating cell means...
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.82 * dense matrix
#> Skipping 'CellMeans': estimated sparse size 1.5 * dense matrix
#> Done!
```
