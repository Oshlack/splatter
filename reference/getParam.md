# Get a parameter

Accessor function for getting parameter values.

## Usage

``` r
getParam(object, name)

# S4 method for class 'Params'
getParam(object, name)
```

## Arguments

- object:

  object to get parameter from.

- name:

  name of the parameter to get.

## Value

The extracted parameter value

## Examples

``` r
params <- newSimpleParams()
getParam(params, "nGenes")
#> [1] 10000
```
