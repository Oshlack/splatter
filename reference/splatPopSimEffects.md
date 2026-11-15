# Add eQTL effects to means matrix

Add eQTL effects and non-eQTL group effects to simulated means matrix.
The eQTL effects are incorporated using the following equation: \$\$Ygs
= (ESg x Mgs x Gs) + Mgs \$\$ Where Ygs is the mean for gene g and
sample s, ESg is the effect size assigned to g, Mgs is the mean
expression assigned to g for s, and Gs is the genotype (number of minor
alleles) for s. Non-eQTL group effects are incorporated as: \$\$Ygs =
Mgs x GEg\$\$ Where GEg is the group effect (i.e. differential
expression) assigned to g. To simulate multiple gene mean matrices with
different group effects, this function can be run with `id` designating
the group id.

## Usage

``` r
splatPopSimEffects(id, key, conditions, vcf, means.pop)
```

## Arguments

- id:

  The group ID (e.g. "global" or "g1")

- key:

  Partial splatPop key data.frame.

- conditions:

  array of condition assignments for each sample

- vcf:

  VariantAnnotation object containing genotypes of samples.

- means.pop:

  Population mean gene expression matrix

## Value

data.frame of gene mean expression levels WITH eQTL effects.
