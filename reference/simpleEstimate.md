# Estimate simple simulation parameters

Estimate simulation parameters for the simple simulation from a real
dataset.

## Usage

``` r
simpleEstimate(counts, params = newSimpleParams())

# S3 method for class 'SingleCellExperiment'
simpleEstimate(counts, params = newSimpleParams())

# S3 method for class 'matrix'
simpleEstimate(counts, params = newSimpleParams())
```

## Arguments

- counts:

  either a counts matrix or a SingleCellExperiment object containing
  count data to estimate parameters from.

- params:

  SimpleParams object to store estimated values in.

## Value

SimpleParams object containing the estimated parameters.

## Details

The `nGenes` and `nCells` parameters are taken from the size of the
input data. The mean parameters are estimated by fitting a gamma
distribution to the library size normalised mean expression level using
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).
See
[`SimpleParams`](http://oshlacklab.com/splatter/reference/SimpleParams.md)
for more details on the parameters.

## Examples

``` r
# Load example data
library(scuttle)
set.seed(1)
sce <- mockSCE()

params <- simpleEstimate(sce)
params
#> A Params object of class SimpleParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>    2000      200   787110 
#> 
#> 3 additional parameters 
#> 
#> Mean: 
#>              (RATE)              (SHAPE) 
#> 0.00315146081285586    0.562634238245172 
#> 
#> Counts: 
#> [Dispersion] 
#>          0.1 
#> 
```
