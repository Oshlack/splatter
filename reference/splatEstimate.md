# Estimate Splat simulation parameters

Estimate simulation parameters for the Splat simulation from a real
dataset. See the individual estimation functions for more details on how
this is done.

## Usage

``` r
splatEstimate(counts, params = newSplatParams())

# S3 method for class 'SingleCellExperiment'
splatEstimate(counts, params = newSplatParams())

# S3 method for class 'matrix'
splatEstimate(counts, params = newSplatParams())
```

## Arguments

- counts:

  either a counts matrix or a SingleCellExperiment object containing
  count data to estimate parameters from.

- params:

  SplatParams object to store estimated values in.

## Value

SplatParams object with estimated values.

## See also

[`splatEstMean`](http://oshlacklab.com/splatter/reference/splatEstMean.md),
[`splatEstLib`](http://oshlacklab.com/splatter/reference/splatEstLib.md),
[`splatEstOutlier`](http://oshlacklab.com/splatter/reference/splatEstOutlier.md),
[`splatEstBCV`](http://oshlacklab.com/splatter/reference/splatEstBCV.md),
[`splatEstDropout`](http://oshlacklab.com/splatter/reference/splatEstDropout.md)

## Examples

``` r
# Load example data
library(scuttle)
set.seed(1)
sce <- mockSCE()

params <- splatEstimate(sce)
#> NOTE: Library sizes have been found to be normally distributed instead of log-normal. You may want to check this is correct.
params
#> A Params object of class SplatParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>    2000      200   787110 
#> 
#> 29 additional parameters 
#> 
#> Batches: 
#>     [BATCHES]  [BATCH CELLS]     [Location]        [Scale]       [Remove] 
#>             1            200            0.1            0.1          FALSE 
#> 
#> Mean: 
#>              (RATE)              (SHAPE) 
#> 0.00296241290380408    0.496975675076379 
#> 
#> Library size: 
#>       (LOCATION)           (SCALE)            (NORM) 
#>       357331.235  11607.2332293176              TRUE 
#> 
#> Exprs outliers: 
#> (PROBABILITY)     (Location)        (Scale) 
#>             0              4            0.5 
#> 
#> Groups: 
#>      [Groups]  [Group Probs] 
#>             1              1 
#> 
#> Diff expr: 
#> [Probability]    [Down Prob]     [Location]        [Scale] 
#>           0.1            0.5            0.1            0.4 
#> 
#> BCV: 
#>     (COMMON DISP)              (DOF) 
#> 0.752043426792849                Inf 
#> 
#> Dropout: 
#>           [Type]        (MIDPOINT)           (SHAPE) 
#>             none  2.71153553670808  -1.3720928186731 
#> 
#> Paths: 
#>         [From]         [Steps]          [Skew]    [Non-linear]  [Sigma Factor] 
#>              0             100             0.5             0.1             0.8 
#> 
```
