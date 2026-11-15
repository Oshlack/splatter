# New Params

Create a new Params object. Functions exist for each of the different
Params subtypes.

## Usage

``` r
newBASiCSParams(...)

newKersplatParams(...)

newLun2Params(...)

newLunParams(...)

newPhenoParams(...)

newSCDDParams(...)

newSimpleParams(...)

newSparseDCParams(...)

newSplatParams(...)

newSplatPopParams(...)

newZINBParams(...)
```

## Arguments

- ...:

  additional parameters passed to
  [`setParams`](http://oshlacklab.com/splatter/reference/setParams.md).

## Value

New Params object.

## Examples

``` r
params <- newSimpleParams()
params <- newSimpleParams(nGenes = 200, nCells = 10)
```
