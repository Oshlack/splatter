# Assign eGenes-eSNPs pairs and effect sizes.

Randomly pairs N genes (eGene) a SNP (eSNP) within the window size
(eqtl.dist) and assigns each pair an effect size sampled from a gamma
distribution parameterized using the effect sizes from a real eQTL
study.

## Usage

``` r
splatPopeQTLEffects(params, key, vcf)
```

## Arguments

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- key:

  Partial splatPop key data.frame.

- vcf:

  VariantAnnotation object containing genotypes of samples.

## Value

The key updated with assigned eQTL effects.
