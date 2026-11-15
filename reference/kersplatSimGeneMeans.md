# Simulate Kersplat gene means

Simulate Kersplat gene means

## Usage

``` r
kersplatSimGeneMeans(params, verbose)
```

## Arguments

- params:

  KersplatParams object containing simulation parameters.

- verbose:

  logical. Whether to print progress messages

## Value

KersplatParams object with gene means

## Details

Gene means are simulated in one of two ways depending on the value of
the `mean.method` parameter.

If `mean.method` is "fit" (default) then means are sampled from a Gamma
distribution with shape equals `mean.shape` and rate equals `mean.rate`.
Expression outliers are then added by replacing some values with the
median multiplied by a factor from a log-normal distribution. This is
the same process used for the Splat simulation.

If `mean.method` is "density" then means are sampled from the density
object in the `mean.density` parameter using a rejection sampling
method. This approach is more flexible but may violate some statistical
assumptions.
