# Get parameters

Get multiple parameter values from a Params object.

## Usage

``` r
getParams(params, names)
```

## Arguments

- params:

  Params object to get values from.

- names:

  vector of names of the parameters to get.

## Value

List with the values of the selected parameters.

## Examples

``` r
params <- newSimpleParams()
getParams(params, c("nGenes", "nCells", "mean.rate"))
#> $nGenes
#> [1] 10000
#> 
#> $nCells
#> [1] 100
#> 
#> $mean.rate
#> [1] 0.3
#> 
```
