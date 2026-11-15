# PhenoPath simulation

Simulate counts from a pseudotime trajectory using the PhenoPath method.

## Usage

``` r
phenoSimulate(params = newPhenoParams(), sparsify = TRUE, verbose = TRUE, ...)
```

## Arguments

- params:

  PhenoParams object containing simulation parameters.

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
[`simulate_phenopath`](https://rdrr.io/pkg/phenopath/man/simulate_phenopath.html)
that takes a
[`PhenoParams`](http://oshlacklab.com/splatter/reference/PhenoParams.md),
runs the simulation then converts the output from log-expression to
counts and returns a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object. The original simulated log-expression values are returned in the
`LogExprs` assay. See
[`simulate_phenopath`](https://rdrr.io/pkg/phenopath/man/simulate_phenopath.html)
and the PhenoPath paper for more details about how the simulation works.

## References

Campbell K, Yau C. Uncovering genomic trajectories with heterogeneous
genetic and environmental backgrounds across single-cells and
populations. bioRxiv (2017).

Paper:
[10.1101/159913](http://oshlacklab.com/splatter/reference/10.1101/159913)

Code: <https://github.com/kieranrcampbell/phenopath>

## Examples

``` r
if (requireNamespace("phenopath", quietly = TRUE)) {
    sim <- phenoSimulate()
}
#> Simulating counts...
#> Warning: `as_data_frame()` was deprecated in tibble 2.0.0.
#> ℹ Please use `as_tibble()` (with slightly different semantics) to convert to a
#>   tibble, or `as.data.frame()` to convert to a data frame.
#> ℹ The deprecated feature was likely used in the phenopath package.
#>   Please report the issue to the authors.
#> Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if
#> `.name_repair` is omitted as of tibble 2.0.0.
#> ℹ Using compatibility `.name_repair`.
#> ℹ The deprecated feature was likely used in the tibble package.
#>   Please report the issue at <https://github.com/tidyverse/tibble/issues>.
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.52 * dense matrix
#> Skipping 'LogExprs': estimated sparse size 1.5 * dense matrix
```
