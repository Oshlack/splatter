# Assign group-specific eQTL and DEGs.

If groups \> 1, n eSNP-eGene pairs (n = 'eqtl.group.specific') are
randomly assigned as group specific.

## Usage

``` r
splatPopGroupEffects(params, key, groups)
```

## Arguments

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- key:

  Partial splatPop key data.frame.

- groups:

  array of group names

## Value

The key updated with group eQTL and DE effects.
