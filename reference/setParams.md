# Set parameters

Set multiple parameters in a Params object.

## Usage

``` r
setParams(object, update = NULL, ...)

# S4 method for class 'KersplatParams'
setParams(object, update = NULL, ...)

# S4 method for class 'Params'
setParams(object, update = NULL, ...)

# S4 method for class 'SplatParams'
setParams(object, update = NULL, ...)
```

## Arguments

- object:

  Params object to set parameters in.

- update:

  list of parameters to set where `names(update)` are the names of the
  parameters to set and the items in the list are values.

- ...:

  additional parameters to set. These are combined with any parameters
  specified in `update`.

## Value

Params object with updated values.

## Details

Each parameter is set by a call to
[`setParam`](http://oshlacklab.com/splatter/reference/setParam.md). If
the same parameter is specified multiple times it will be set multiple
times. Parameters can be specified using a list via `update` (useful
when collecting parameter values in some way) or individually (useful
when setting them manually), see examples.

## Examples

``` r
params <- newSimpleParams()
params
#> A Params object of class SimpleParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (Genes)  (Cells)   [SEED] 
#>   10000      100   425263 
#> 
#> 3 additional parameters 
#> 
#> Mean: 
#>  (Rate)  (Shape) 
#>     0.3      0.4 
#> 
#> Counts: 
#> [Dispersion] 
#>          0.1 
#> 
# Set individually
params <- setParams(params, nGenes = 1000, nCells = 50)
params
#> A Params object of class SimpleParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>    1000       50   425263 
#> 
#> 3 additional parameters 
#> 
#> Mean: 
#>  (Rate)  (Shape) 
#>     0.3      0.4 
#> 
#> Counts: 
#> [Dispersion] 
#>          0.1 
#> 
# Set via update list
params <- setParams(params, list(mean.rate = 0.2, mean.shape = 0.8))
params
#> A Params object of class SimpleParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>    1000       50   425263 
#> 
#> 3 additional parameters 
#> 
#> Mean: 
#>  (RATE)  (SHAPE) 
#>     0.2      0.8 
#> 
#> Counts: 
#> [Dispersion] 
#>          0.1 
#> 
```
