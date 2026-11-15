# Expand parameters

Expand the parameters that can be vectors so that they are the same
length as the number of groups. Work is done by `paramsExpander` called
from each method. Expansions are stored using
[`setParamsUnchecked`](http://oshlacklab.com/splatter/reference/setParamsUnchecked.md).

## Usage

``` r
expandParams(object, ...)

# S4 method for class 'BASiCSParams'
expandParams(object)

# S4 method for class 'LunParams'
expandParams(object)

# S4 method for class 'Params'
expandParams(object, vectors, n)

# S4 method for class 'SplatParams'
expandParams(object)

# S4 method for class 'SplatPopParams'
expandParams(object)

paramsExpander(object, vectors, n)
```

## Arguments

- object:

  object to expand.

- ...:

  additional arguments.

- vectors:

  names of vector parameters to expand

- n:

  number of times to repeat each parameter

## Value

Expanded object.
