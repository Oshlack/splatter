# Simulate Kersplat final counts

Simulate the final counts matrix for a Kersplat simulation

## Usage

``` r
kersplatSimCounts(sim, params, verbose)
```

## Arguments

- sim:

  SingleCellExperiment containing simulation.

- params:

  KersplatParams object with simulation parameters.

- verbose:

  logical. Whether to print progress messages

## Value

SingleCellExperiment with counts matrix

## Details

The cell counts matrix and ambient counts matrix are added together. The
result is then downsampled to the cell library size (for cells and
doublets) or the ambient library size (for empty cells) using the
[`downsampleMatrix`](https://rdrr.io/pkg/scuttle/man/downsampleMatrix.html)
function.

## See also

[`downsampleMatrix`](https://rdrr.io/pkg/scuttle/man/downsampleMatrix.html)
