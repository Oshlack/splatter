# Estimate Lun simulation parameters

Estimate simulation parameters for the Lun simulation from a real
dataset.

## Usage

``` r
lunEstimate(counts, params = newLunParams())

# S3 method for class 'SingleCellExperiment'
lunEstimate(counts, params = newLunParams())

# S3 method for class 'matrix'
lunEstimate(counts, params = newLunParams())
```

## Arguments

- counts:

  either a counts matrix or an SingleCellExperiment object containing
  count data to estimate parameters from.

- params:

  LunParams object to store estimated values in.

## Value

LunParams object containing the estimated parameters.

## Details

The `nGenes` and `nCells` parameters are taken from the size of the
input data. No other parameters are estimated. See
[`LunParams`](http://oshlacklab.com/splatter/reference/LunParams.md) for
more details on the parameters.

## Examples

``` r
# Load example data
library(scuttle)
set.seed(1)
sce <- mockSCE()

params <- lunEstimate(sce)
params
#> A Params object of class LunParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>    2000      200   787110 
#> 
#> 9 additional parameters 
#> 
#> Mean: 
#>  [Rate]  [Shape] 
#>       2        2 
#> 
#> Counts: 
#> [Dispersion] 
#>          0.1 
#> 
#> Groups: 
#>      [GROUPS]  [GROUP CELLS] 
#>             1            200 
#> 
#> Diff expr: 
#>   [Genes]  [Up Prop]    [Up FC]  [Down FC] 
#>      1000        0.5          5          0 
#> 
```
