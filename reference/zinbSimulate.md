# ZINB-WaVE simulation

Simulate counts using the ZINB-WaVE method.

## Usage

``` r
zinbSimulate(params = newZINBParams(), sparsify = TRUE, verbose = TRUE, ...)
```

## Arguments

- params:

  ZINBParams object containing simulation parameters.

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
[`zinbSim`](https://rdrr.io/pkg/zinbwave/man/zinbSim.html) that takes a
[`ZINBParams`](http://oshlacklab.com/splatter/reference/ZINBParams.md),
runs the simulation then converts the output to a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object. See [`zinbSim`](https://rdrr.io/pkg/zinbwave/man/zinbSim.html)
and the ZINB-WaVE paper for more details about how the simulation works.

## References

Campbell K, Yau C. Uncovering genomic trajectories with heterogeneous
genetic and environmental backgrounds across single-cells and
populations. bioRxiv (2017).

Risso D, Perraudeau F, Gribkova S, Dudoit S, Vert J-P. ZINB-WaVE: A
general and flexible method for signal extraction from single-cell
RNA-seq data bioRxiv (2017).

Paper:
[10.1101/125112](http://oshlacklab.com/splatter/reference/10.1101/125112)

Code: <https://github.com/drisso/zinbwave>

## Examples

``` r
if (requireNamespace("zinbwave", quietly = TRUE)) {
    sim <- zinbSimulate()
}
#> Simulating counts...
#> Creating final dataset...
#> Sparsifying assays...
#> Automatically converting to sparse matrices, threshold = 0.95
#> Converting 'counts' to sparse matrix: estimated sparse size 0.37 * dense matrix
#> Converting 'TrueCounts' to sparse matrix: estimated sparse size 0.74 * dense matrix
#> Skipping 'Dropouts': estimated sparse size 1.54 * dense matrix
```
