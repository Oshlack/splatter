# Kersplat setup

Setup the parameters required for the Kersplat simulation

## Usage

``` r
kersplatSetup(params = newKersplatParams(), verbose = TRUE, ...)
```

## Arguments

- params:

  KersplatParams object containing simulation parameters.

- verbose:

  logical. Whether to print progress messages

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

A complete KersplatParams object

## Details

The first stage is a two-step Kersplat simulation is to generate some of
the intermediate parameters. The resulting parameters allow multiple
simulated datasets to be generated from the same biological structure
(using
[`kersplatSample`](http://oshlacklab.com/splatter/reference/kersplatSample.md)).
As with all the other parameters these values can be manually
overwritten if desired.

The setup involves the following steps:

1.  Generate a gene network (if not already present)

2.  Select regulator genes (if not already present)

3.  Simulate gene means (if not already present)

4.  Simulate cell paths

The resulting
[`KersplatParams`](http://oshlacklab.com/splatter/reference/KersplatParams.md)
object will have the following parameters set (if they weren't already).

- `mean.values`

- `network.graph`

- `network.regsSet`

- `paths.means`

See
[`KersplatParams`](http://oshlacklab.com/splatter/reference/KersplatParams.md)
for more details about these parameters and the functions for the
individual steps for more details about the process.

## See also

[`kersplatGenNetwork`](http://oshlacklab.com/splatter/reference/kersplatGenNetwork.md),
[`kersplatSelectRegs`](http://oshlacklab.com/splatter/reference/kersplatSelectRegs.md),
[`kersplatSimGeneMeans`](http://oshlacklab.com/splatter/reference/kersplatSimGeneMeans.md),
[`kersplatSimPaths`](http://oshlacklab.com/splatter/reference/kersplatSimPaths.md),
[`KersplatParams`](http://oshlacklab.com/splatter/reference/KersplatParams.md)

## Examples

``` r

if (requireNamespace("igraph", quietly = TRUE)) {
    params <- kersplatSetup()
}
#> Setting up parameters...
#> Generating gene network...
#> Selecting regulators...
#> Simulating means...
#> Sampling from gamma distribution...
#> Simulating paths...
#> Simulating path 1...
```
