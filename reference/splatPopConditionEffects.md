# Assign Condition-specific eQTL and DEGs.

If nConditions \> 1, n eSNP-eGene pairs (n = 'eqtl.condition.specific')
are randomly assigned as condition specific.

## Usage

``` r
splatPopConditionEffects(params, key, conditions)
```

## Arguments

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- key:

  Partial splatPop key data.frame.

- conditions:

  array of condition names

## Value

The key updated with conditional eQTL and DE effects.
