# Set a parameter

Function for setting parameter values.

## Usage

``` r
setParam(object, name, value)

# S4 method for class 'BASiCSParams'
setParam(object, name, value)

# S4 method for class 'KersplatParams'
setParam(object, name, value)

# S4 method for class 'Lun2Params'
setParam(object, name, value)

# S4 method for class 'LunParams'
setParam(object, name, value)

# S4 method for class 'Params'
setParam(object, name, value)

# S4 method for class 'PhenoParams'
setParam(object, name, value)

# S4 method for class 'SCDDParams'
setParam(object, name, value)

# S4 method for class 'SplatParams'
setParam(object, name, value)

# S4 method for class 'SplatPopParams'
setParam(object, name, value)

# S4 method for class 'ZINBParams'
setParam(object, name, value)
```

## Arguments

- object:

  object to set parameter in.

- name:

  name of the parameter to set.

- value:

  value to set the parameter to.

## Value

Object with new parameter value.

## Examples

``` r
params <- newSimpleParams()
setParam(params, "nGenes", 100)
#> A Params object of class SimpleParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (Cells)   [SEED] 
#>     100      100   316949 
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
```
