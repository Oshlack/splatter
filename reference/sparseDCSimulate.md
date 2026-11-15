# SparseDC simulation

Simulate counts from cluster in two conditions using the SparseDC
method.

## Usage

``` r
sparseDCSimulate(
  params = newSparseDCParams(),
  sparsify = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- params:

  SparseDCParams object containing simulation parameters.

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
[`sim_data`](https://rdrr.io/pkg/SparseDC/man/sim_data.html) that takes
a
[`SparseDCParams`](http://oshlacklab.com/splatter/reference/SparseDCParams.md),
runs the simulation then converts the output from log-expression to
counts and returns a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object. The original simulated log-expression values are returned in the
`LogExprs` assay. See
[`sim_data`](https://rdrr.io/pkg/SparseDC/man/sim_data.html) and the
SparseDC paper for more details about how the simulation works.

## References

Campbell K, Yau C. Uncovering genomic trajectories with heterogeneous
genetic and environmental backgrounds across single-cells and
populations. bioRxiv (2017).

Barron M, Zhang S, Li J. A sparse differential clustering algorithm for
tracing cell type changes via single-cell RNA-sequencing data. Nucleic
Acids Research (2017).

Paper:
[10.1093/nar/gkx1113](http://oshlacklab.com/splatter/reference/10.1093/nar/gkx1113)

## Examples

``` r
if (requireNamespace("SparseDC", quietly = TRUE)) {
    sim <- sparseDCSimulate()
}
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.58 * dense matrix
#> Skipping 'LogExprs': estimated sparse size 1.5 * dense matrix
```
