# Simulate library sizes

Simulate expected library sizes. Typically a log-normal distribution is
used but there is also the option to use a normal distribution. In this
case any negative values are set to half the minimum non-zero value.

## Usage

``` r
splatSimLibSizes(sim, params)
```

## Arguments

- sim:

  SingleCellExperiment to add library size to.

- params:

  SplatParams object with simulation parameters.

## Value

SingleCellExperiment with simulated library sizes.
