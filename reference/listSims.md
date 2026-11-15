# List simulations

List all the simulations that are currently available in Splatter with a
brief description.

## Usage

``` r
listSims(print = TRUE)
```

## Arguments

- print:

  logical. Whether to print to the console.

## Value

Invisibly returns a data.frame containing the information that is
displayed.

## Examples

``` r
listSims()
#> Splatter currently contains 14 simulations 
#> 
#> Splat (splat) 
#> DOI: 10.1186/s13059-017-1305-0   GitHub: Oshlack/splatter   Dependencies:  
#> The Splat simulation generates means from a gamma distribution, adjusts them for BCV and generates counts from a gamma-poisson. Dropout and batch effects can be optionally added. 
#> 
#> Splat Single (splatSingle) 
#> DOI: 10.1186/s13059-017-1305-0   GitHub: Oshlack/splatter   Dependencies:  
#> The Splat simulation with a single population. 
#> 
#> Splat Groups (splatGroups) 
#> DOI: 10.1186/s13059-017-1305-0   GitHub: Oshlack/splatter   Dependencies:  
#> The Splat simulation with multiple groups. Each group can have it's own differential expression probability and fold change distribution. 
#> 
#> Splat Paths (splatPaths) 
#> DOI: 10.1186/s13059-017-1305-0   GitHub: Oshlack/splatter   Dependencies:  
#> The Splat simulation with differentiation paths. Each path can have it's own length, skew and probability. Genes can change in non-linear ways. 
#> 
#> Kersplat (kersplat) 
#> DOI:      GitHub: Oshlack/splatter   Dependencies: scuttle, igraph 
#> The Kersplat simulation extends the Splat model by adding a gene network, more complex cell structure, doublets and empty cells (Experimental). 
#> 
#> splatPop (splatPop) 
#> DOI: 10.1186/s13059-021-02546-1      GitHub: Oshlack/splatter   Dependencies: VariantAnnotation, preprocessCore 
#> The splatPop simulation enables splat simulations to be generated for multiple individuals in a population, accounting for correlation structure by simulating expression quantitative trait loci (eQTL). 
#> 
#> Simple (simple) 
#> DOI: 10.1186/s13059-017-1305-0   GitHub: Oshlack/splatter   Dependencies:  
#> A simple simulation with gamma means and negative binomial counts. 
#> 
#> Lun (lun) 
#> DOI: 10.1186/s13059-016-0947-7   GitHub: MarioniLab/Deconvolution2016   Dependencies:  
#> Gamma distributed means and negative binomial counts. Cells are given a size factor and differential expression can be simulated with fixed fold changes. 
#> 
#> Lun 2 (lun2) 
#> DOI: 10.1093/biostatistics/kxw055    GitHub: MarioniLab/PlateEffects2016    Dependencies: scran, scuttle, lme4, pscl, limSolve 
#> Negative binomial counts where the means and dispersions have been sampled from a real dataset. The core feature of the Lun 2 simulation is the addition of plate effects. Differential expression can be added between two groups of plates and optionally a zero-inflated negative-binomial can be used. 
#> 
#> scDD (scDD) 
#> DOI: 10.1186/s13059-016-1077-y   GitHub: kdkorthauer/scDD   Dependencies: scDD 
#> The scDD simulation samples a given dataset and can simulate differentially expressed and differentially distributed genes between two conditions. 
#> 
#> BASiCS (BASiCS) 
#> DOI: 10.1371/journal.pcbi.1004333    GitHub: catavallejos/BASiCS    Dependencies: BASiCS 
#> The BASiCS simulation is based on a bayesian model used to deconvolve biological and technical variation and includes spike-ins and batch effects. 
#> 
#> PhenoPath (pheno) 
#> DOI: 10.1038/s41467-018-04696-6      GitHub: kieranrcampbell/phenopath      Dependencies: phenopath 
#> The PhenoPath simulation produces a pseudotime trajectory with different types of genes. 
#> 
#> ZINB-WaVE (zinb) 
#> DOI: 10.1038/s41467-017-02554-5      GitHub: drisso/zinbwave    Dependencies: zinbwave 
#> The ZINB-WaVE simulation simulates counts from a sophisticated zero-inflated negative-binomial distribution including cell and gene-level covariates. 
#> 
#> SparseDC (sparseDC) 
#> DOI: 10.1093/nar/gkx1113     GitHub: cran/SparseDC      Dependencies: SparseDC 
#> The SparseDC simulation simulates a set of clusters across two conditions, where some clusters may be present in only one condition. 
#> 
```
