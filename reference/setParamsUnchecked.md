# Set parameters UNCHECKED

Set multiple parameters in a Params object.

## Usage

``` r
setParamsUnchecked(params, update = NULL, ...)
```

## Arguments

- params:

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
when setting them manually), see examples. THE FINAL OBJECT IS NOT
CHECKED FOR VALIDITY!
