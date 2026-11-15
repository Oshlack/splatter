# Simulate true counts

Simulate a true counts matrix. Counts are simulated from a poisson
distribution where Each gene in each cell has it's own mean based on the
group (or path position), expected library size and BCV.

## Usage

``` r
splatSimTrueCounts(sim, params)
```

## Arguments

- sim:

  SingleCellExperiment to add true counts to.

- params:

  SplatParams object with simulation parameters.

## Value

SingleCellExperiment with simulated true counts.
