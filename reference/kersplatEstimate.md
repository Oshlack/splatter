# Estimate Kersplat simulation parameters

Estimate simulation parameters for the Kersplat simulation from a real
dataset. See the individual estimation functions for more details on how
this is done.

## Usage

``` r
kersplatEstimate(counts, params = newKersplatParams(), verbose = TRUE)

# S3 method for class 'SingleCellExperiment'
kersplatEstimate(counts, params = newKersplatParams(), verbose = TRUE)

# S3 method for class 'matrix'
kersplatEstimate(counts, params = newKersplatParams(), verbose = TRUE)
```

## Arguments

- counts:

  either a counts matrix or a SingleCellExperiment object containing
  count data to estimate parameters from.

- params:

  KersplatParams object to store estimated values in.

- verbose:

  logical. Whether to print progress messages.

## Value

KersplatParams object containing the estimated parameters.

## See also

[`kersplatEstMean`](http://oshlacklab.com/splatter/reference/kersplatEstMean.md),
[`kersplatEstBCV`](http://oshlacklab.com/splatter/reference/kersplatEstBCV.md),
[`kersplatEstLib`](http://oshlacklab.com/splatter/reference/kersplatEstLib.md)

## Examples

``` r

if (requireNamespace("igraph", quietly = TRUE)) {
    # Load example data
    library(scuttle)
    set.seed(1)
    sce <- mockSCE()

    params <- kersplatEstimate(sce)
    params
}
#> Warning: The Kersplat simulation is still experimental and may produce unreliable results. Please try it and report any issues to https://github.com/Oshlack/splatter/issues. The development version may have improved features.
#> Estimating mean parameters...
#> Selected MGE (CvM) fit
#> Estimating expression outlier parameters...
#> Estimating BCV parameters...
#> Raw: 2.6081737071714 A: 8.12109522965305 B: 1.31905398252679 C: -8.46037793603958 Y: 1.36293828975133
#> Estimating library size parameters...
#> Selected MGE (CvM) fit
#> A Params object of class KersplatParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>    2000      200   787110 
#> 
#> 24 additional parameters 
#> 
#> Mean: 
#>                  (RATE)                  (SHAPE)               (OUT PROB) 
#>     0.00272353361858179        0.475099389599897                        0 
#>          (Out Location)              (Out Scale)                (DENSITY) 
#>                       4                      0.5  Object of class density 
#>                [Method]                 [Values] 
#>                     fit                  Not set 
#> 
#> BCV: 
#>    (COMMON DISP)             [DoF] 
#> 1.11818825375479                60 
#> 
#> Network: 
#>    [Graph]     [nRegs]   [regsSet] 
#>    Not set         100       FALSE 
#> 
#> Paths: 
#> [nPrograms]      [Means] 
#>          10      Not set 
#> 
#> [Design]
#> data.frame (1 x 3) with columns: Path, From, Steps 
#>   Path From Steps
#> 1    1    0   100
#> 
#> Library size: 
#>              (LOCATION)                  (SCALE)                (DENSITY) 
#>        12.7854772402738       0.0327578186101294  Object of class density 
#>                [Method] 
#>                     fit 
#> 
#> Cells: 
#> 
#> [Design]
#> data.frame (1 x 4) with columns: Path, Probability, Alpha, Beta 
#>   Path Probability Alpha Beta
#> 1    1           1     1    0
#> 
#> Doublets: 
#> [Prop] 
#>      0 
#> 
#> Ambient: 
#> [Scale]  [Empty] 
#>    0.05        0 
#> 
```
