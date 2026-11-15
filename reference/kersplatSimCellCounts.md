# Simulate Kersplat cell counts

Simulate cell counts for the Kersplat simulation

## Usage

``` r
kersplatSimCellCounts(sim, params, verbose)
```

## Arguments

- sim:

  SingleCellExperiment containing simulation.

- params:

  KersplatParams object with simulation parameters.

- verbose:

  logical. Whether to print progress messages

## Value

SingleCellExperiment with cell counts

## Details

Counts are sampled from a Poisson distribution with lambda equal to the
cell means matrix.
