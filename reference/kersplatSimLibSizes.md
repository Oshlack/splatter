# Simulate Kersplat library sizes

Generate library sizes for cells in the Kersplat simulation

## Usage

``` r
kersplatSimLibSizes(sim, params, verbose)
```

## Arguments

- sim:

  SingleCellExperiment containing simulation.

- params:

  KersplatParams object with simulation parameters.

- verbose:

  logical. Whether to print progress messages

## Value

SingleCellExperiment with library sizes

## Details

Library sizes are simulated in one of two ways depending on the value of
the `lib.method` parameter.

If `lib.method` is "fit" (default) then means are sampled from a
log-normal distribution with meanlog equals `lib.loc` and sdlog equals
`lib.scale`.

If `mean.method` is "density" then library sizes are sampled from the
density object in the `lib.density` parameter using a rejection sampling
method. This approach is more flexible but may violate some statistical
assumptions.

Ambient library sizes are also generated from a log-normal distribution
based on the parameters for the cell library size and adjusted using the
`ambient.scale` parameter.
