# Simulate Kersplat ambient counts

Simulate Kersplat ambient counts

## Usage

``` r
kersplatSimAmbientCounts(sim, params, verbose)
```

## Arguments

- sim:

  SingleCellExperiment containing simulation.

- params:

  KersplatParams object with simulation parameters.

- verbose:

  logical. Whether to print progress messages

## Value

SingleCellExperiment with ambient counts

## Details

The overall expression profile to calculated by averaging the cell
counts of the (non-empty) cells. This is then multiplied by the ambient
library sizes to get a mean for each cell. Counts are then sampled from
a Poisson distribution using these means.
