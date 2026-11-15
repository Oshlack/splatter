# Estimate Lun2 simulation parameters

Estimate simulation parameters for the Lun2 simulation from a real
dataset.

## Usage

``` r
lun2Estimate(
  counts,
  plates,
  params = newLun2Params(),
  min.size = 200,
  verbose = TRUE,
  BPPARAM = SerialParam()
)

# S3 method for class 'SingleCellExperiment'
lun2Estimate(
  counts,
  plates,
  params = newLun2Params(),
  min.size = 200,
  verbose = TRUE,
  BPPARAM = SerialParam()
)

# S3 method for class 'matrix'
lun2Estimate(
  counts,
  plates,
  params = newLun2Params(),
  min.size = 200,
  verbose = TRUE,
  BPPARAM = SerialParam()
)
```

## Arguments

- counts:

  either a counts matrix or a SingleCellExperiment object containing
  count data to estimate parameters from.

- plates:

  integer vector giving the plate that each cell originated from.

- params:

  Lun2Params object to store estimated values in.

- min.size:

  minimum size of clusters when identifying group of cells in the data.

- verbose:

  logical. Whether to show progress messages.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance giving the parallel back-end to be used. Default is
  [`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  which uses a single core.

## Value

LunParams object containing the estimated parameters.

## Details

See
[`Lun2Params`](http://oshlacklab.com/splatter/reference/Lun2Params.md)
for more details on the parameters.

## Examples

``` r
# \donttest{
# Load example data
library(scuttle)
set.seed(1)
sce <- mockSCE()

plates <- as.numeric(factor(colData(sce)$Mutation_Status))
params <- lun2Estimate(sce, plates, min.size = 20)
#> Estimating number of groups...
#> Computing normalisation factors...
#> Estimating dispersions...
#> Estimating gene means...
#> Estimating plate effects...
#> Estimating zero-inflated parameters...
#> Warning: NaNs produced
#> Warning: NaNs produced
params
#> A Params object of class Lun2Params 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>    2000      200   787110 
#> 
#> 11 additional parameters 
#> 
#> Genes: 
#> 
#> (PARAMS)
#> data.frame (2000 x 2) with columns: Mean, Disp 
#>         Mean      Disp
#> 1  61.200373  2.391091
#> 2   3.915675 19.460132
#> 3 218.990476  0.998549
#> 4   6.268318 25.938766
#> # ... with 1996 more rows
#> 
#> (ZI PARAMS)
#> data.frame (2000 x 3) with columns: Mean, Disp, Prop 
#>         Mean      Disp         Prop
#> 1  26.795530  4.629399 3.119853e-06
#> 2   2.376801 19.503446 7.215299e-05
#> 3 218.990476  0.998549 0.000000e+00
#> 4   3.194146 25.610243 4.154911e-04
#> # ... with 1996 more rows
#> 
#> Plates: 
#>         [NUMBER]        [Modifier]        (VARIANCE) 
#>                2                 1  19.8392380729304 
#> 
#> Cells: 
#>                           [PLATES]                     (LIBRARY SIZES) 
#>                     1, 2, 1, 1,...  344015, 334664, 367199, 349098,... 
#>                     [Lib Size Mod] 
#>                                  1 
#> 
#> Diff Expr: 
#>       [Genes]  [Fold change] 
#>             0              3 
#> 
# }
```
