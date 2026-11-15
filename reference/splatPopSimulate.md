# splatPop simulation

Simulate scRNA-seq count data using the splat model for a population of
individuals with correlation structure.

## Usage

``` r
splatPopSimulate(
  params = newSplatPopParams(nGenes = 50),
  vcf = mockVCF(),
  method = c("single", "groups", "paths"),
  gff = NULL,
  eqtl = NULL,
  means = NULL,
  key = NULL,
  counts.only = FALSE,
  sparsify = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- params:

  SplatPopParams object containing parameters for population scale
  simulations. See
  [`SplatPopParams`](http://oshlacklab.com/splatter/reference/SplatPopParams.md)
  for details.

- vcf:

  VariantAnnotation object containing genotypes of samples.

- method:

  which simulation method to use. Options are "single" which produces a
  single population, "groups" which produces distinct groups (eg. cell
  types), "paths" which selects cells from continuous trajectories (eg.
  differentiation processes).

- gff:

  Either NULL or a data.frame object containing a GFF/GTF file.

- eqtl:

  Either NULL or if simulating population parameters directly from
  empirical data, a data.frame with empirical/desired eQTL results. To
  see required format, run
  [`mockEmpiricalSet()`](http://oshlacklab.com/splatter/reference/mockEmpiricalSet.md)
  and see eqtl output.

- means:

  Either NULL or if simulating population parameters directly from
  empirical data, a Matrix of real gene means across a population, where
  each row is a gene and each column is an individual in the population.
  To see required format, run
  [`mockEmpiricalSet()`](http://oshlacklab.com/splatter/reference/mockEmpiricalSet.md)
  and see means output.

- key:

  Either NULL or a data.frame object containing a full or partial
  splatPop key.

- counts.only:

  logical. Whether to save only counts in sce object.

- sparsify:

  logical. Whether to automatically convert assays to sparse matrices if
  there will be a size reduction.

- verbose:

  logical. Whether to print progress messages.

- ...:

  any additional parameter settings to override what is provided in
  `params`.

## Value

SingleCellExperiment object containing simulated counts, intermediate
values like the gene means simulated in `splatPopSimulateMeans`, and
information about the differential expression and eQTL effects assigned
to each gene.

## Details

This functions is for simulating data in a single step. It consists of a
call to
[`splatPopSimulateMeans`](http://oshlacklab.com/splatter/reference/splatPopSimulateMeans.md),
which simulates a mean expression level per gene per sample, followed by
a call to
[`splatPopSimulateSC`](http://oshlacklab.com/splatter/reference/splatPopSimulateSC.md),
which uses the splat model to simulate single-cell counts per
individual. Please see the documentation for those functions for more
details.

## See also

[`splatPopSimulateMeans`](http://oshlacklab.com/splatter/reference/splatPopSimulateMeans.md),
[`splatPopSimulateSC`](http://oshlacklab.com/splatter/reference/splatPopSimulateSC.md)

## Examples

``` r
# \donttest{
if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
    requireNamespace("preprocessCore", quietly = TRUE)) {
    vcf <- mockVCF()
    gff <- mockGFF()
    sim <- splatPopSimulate(vcf = vcf, gff = gff, sparsify = FALSE)
}
#> Designing population...
#> Simulating data for genes in GFF...
#> Simulating gene means for population...
#> Simulating sc counts for Group1...
#> Done!
# }
```
