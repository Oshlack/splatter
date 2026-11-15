# Add conditional DE effects to means matrix

Add conditional DE effects to means matrix

## Usage

``` r
splatPopConditionalEffects(id, key, vcf, means.pop)
```

## Arguments

- id:

  The group ID (e.g. "global" or "g1")

- key:

  Partial splatPop key data.frame.

- vcf:

  VariantAnnotation object containing genotypes of samples.

- means.pop:

  Population mean gene expression matrix

## Value

data.frame of gene mean expression levels WITH eQTL effects.
