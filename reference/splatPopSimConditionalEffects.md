# Add conditional DE effects to means matrix

Add conditional DE effects to means matrix

## Usage

``` r
splatPopSimConditionalEffects(key, means.pop, conditions)
```

## Arguments

- key:

  Partial splatPop key data.frame.

- means.pop:

  matrix or list of matrices with gene means.

- conditions:

  array of condition assignments for each sample

## Value

data.frame of gene mean expression levels WITH conditional DE effects.
