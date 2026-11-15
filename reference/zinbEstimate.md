# Estimate ZINB-WaVE simulation parameters

Estimate simulation parameters for the ZINB-WaVE simulation from a real
dataset.

## Usage

``` r
zinbEstimate(
  counts,
  design.samples = NULL,
  design.genes = NULL,
  common.disp = TRUE,
  iter.init = 2,
  iter.opt = 25,
  stop.opt = 1e-04,
  params = newZINBParams(),
  verbose = TRUE,
  BPPARAM = SerialParam(),
  ...
)

# S3 method for class 'SingleCellExperiment'
zinbEstimate(
  counts,
  design.samples = NULL,
  design.genes = NULL,
  common.disp = TRUE,
  iter.init = 2,
  iter.opt = 25,
  stop.opt = 1e-04,
  params = newZINBParams(),
  verbose = TRUE,
  BPPARAM = SerialParam(),
  ...
)

# S3 method for class 'matrix'
zinbEstimate(
  counts,
  design.samples = NULL,
  design.genes = NULL,
  common.disp = TRUE,
  iter.init = 2,
  iter.opt = 25,
  stop.opt = 1e-04,
  params = newZINBParams(),
  verbose = TRUE,
  BPPARAM = SerialParam(),
  ...
)
```

## Arguments

- counts:

  either a counts matrix or a SingleCellExperiment object containing
  count data to estimate parameters from.

- design.samples:

  design matrix of sample-level covariates.

- design.genes:

  design matrix of gene-level covariates.

- common.disp:

  logical. Whether or not a single dispersion for all features is
  estimated.

- iter.init:

  number of iterations to use for initialization.

- iter.opt:

  number of iterations to use for optimization.

- stop.opt:

  stopping criterion for optimization.

- params:

  ZINBParams object to store estimated values in.

- verbose:

  logical. Whether to print progress messages.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance giving the parallel back-end to be used. Default is
  [`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  which uses a single core.

- ...:

  additional arguments passes to
  [`zinbFit`](https://rdrr.io/pkg/zinbwave/man/zinbFit.html).

## Value

ZINBParams object containing the estimated parameters.

## Details

The function is a wrapper around
[`zinbFit`](https://rdrr.io/pkg/zinbwave/man/zinbFit.html) that takes
the fitted model and inserts it into a
[`ZINBParams`](http://oshlacklab.com/splatter/reference/ZINBParams.md)
object. See
[`ZINBParams`](http://oshlacklab.com/splatter/reference/ZINBParams.md)
for more details on the parameters and
[`zinbFit`](https://rdrr.io/pkg/zinbwave/man/zinbFit.html) for details
of the estimation procedure.

## Examples

``` r
if (requireNamespace("zinbwave", quietly = TRUE)) {
    library(scuttle)
    set.seed(1)
    sce <- mockSCE(ncells = 20, ngenes = 100)

    params <- zinbEstimate(sce)
    params
}
#> Removing all zero genes...
#> Fitting model...
#> Create model:
#> ok
#> Initialize parameters:
#> ok
#> Optimize parameters:
#> Iteration 1
#> penalized log-likelihood = -9797.17477826626
#> After dispersion optimization = -9593.51952142957
#>    user  system elapsed 
#>   0.168   0.000   0.168 
#> After right optimization = -9289.80805409958
#> After orthogonalization = -9289.80805409958
#>    user  system elapsed 
#>   0.032   0.000   0.032 
#> After left optimization = -9285.37113359872
#> After orthogonalization = -9285.37113359872
#> Iteration 2
#> penalized log-likelihood = -9285.37113359872
#> After dispersion optimization = -9235.11877812819
#>    user  system elapsed 
#>   0.101   0.000   0.101 
#> After right optimization = -9229.79020324124
#> After orthogonalization = -9229.79020324124
#>    user  system elapsed 
#>   0.030   0.000   0.029 
#> After left optimization = -9229.25017472597
#> After orthogonalization = -9229.25017472597
#> Iteration 3
#> penalized log-likelihood = -9229.25017472597
#> After dispersion optimization = -9228.65439930358
#>    user  system elapsed 
#>   0.089   0.001   0.090 
#> After right optimization = -9228.44379199689
#> After orthogonalization = -9228.44379199689
#>    user  system elapsed 
#>   0.027   0.000   0.028 
#> After left optimization = -9228.40711392283
#> After orthogonalization = -9228.40711392283
#> Iteration 4
#> penalized log-likelihood = -9228.40711392283
#> ok
#> A Params object of class ZINBParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>     100       20   987270 
#> 
#> 1 additional parameters 
#> 
#> Model: 
#> ZinbModel with 100 features, 20 samples, 0 latent factors and 241 parameters 
#> 
#> Parameters of the ZinbModel 
#> 
#> Model Design: 
#>      (SAMPLES)         (Genes) 
#> 1, 1, 1, 1,...  1, 1, 1, 1,... 
#> 
#> Model Offsets: 
#>            (MU)             (PI) 
#> 20 x 100 matrix  20 x 100 matrix 
#> 
#> Model Indices: 
#> (Sample Mu)    (Gene Mu)  (Sample Pi)    (Gene Pi) 
#>           1            1            1            1 
#> 
#> Model Intercepts: 
#> (Sample Mu)    (Gene Mu)  (Sample Pi)    (Gene Pi) 
#>        TRUE         TRUE         TRUE         TRUE 
#> 
#> Model Latent factors: 
#>           (W) 
#> 20 x 0 matrix 
#> 
#> Model Coefficients: 
#>                                                                     (SAMPLE MU)                                                                        (GENE MU)                                                                      (Latent Mu)                                                                      (SAMPLE PI)                                                                        (GENE PI)                                                                      (Latent Pi) 
#> -0.459768351872422, 1.6211941729215, -0.223102891231897, 0.0429479650554433,...       4.26299873076136, 4.29178876307093, 4.44953735745685, 3.86260700475711,...                                                                   0 x 100 matrix   1.90227701503618, -6.94250263635572, -1.07944945541051, -0.430739469052651,...  -0.983739822631623, -1.35075676098436, -1.23285448546658, -1.26060101467272,...                                                                   0 x 100 matrix 
#> 
#> Model Regularisation: 
#>     (Sample Mu)        (GENE MU)      (Sample Pi)        (GENE PI) 
#>               1                5                1                5 
#>        (LATENT)  (Latent coeffs)           (Zeta)          (Logit) 
#>               5                1              100            0.001 
#> 
```
