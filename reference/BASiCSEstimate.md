# Estimate BASiCS simulation parameters

Estimate simulation parameters for the BASiCS simulation from a real
dataset.

## Usage

``` r
BASiCSEstimate(
  counts,
  spike.info = NULL,
  batch = NULL,
  n = 20000,
  thin = 10,
  burn = 5000,
  regression = TRUE,
  params = newBASiCSParams(),
  verbose = TRUE,
  progress = TRUE,
  ...
)

# S3 method for class 'SingleCellExperiment'
BASiCSEstimate(
  counts,
  spike.info = NULL,
  batch = NULL,
  n = 20000,
  thin = 10,
  burn = 5000,
  regression = TRUE,
  params = newBASiCSParams(),
  verbose = TRUE,
  progress = TRUE,
  ...
)

# S3 method for class 'matrix'
BASiCSEstimate(
  counts,
  spike.info = NULL,
  batch = NULL,
  n = 20000,
  thin = 10,
  burn = 5000,
  regression = TRUE,
  params = newBASiCSParams(),
  verbose = TRUE,
  progress = TRUE,
  ...
)
```

## Arguments

- counts:

  either a counts matrix or a SingleCellExperiment object containing
  count data to estimate parameters from.

- spike.info:

  data.frame describing spike-ins with two columns: "Name" giving the
  names of the spike-in features (must match `rownames(counts)`) and
  "Input" giving the number of input molecules.

- batch:

  vector giving the batch that each cell belongs to.

- n:

  total number of MCMC iterations. Must be `>= max(4, thin)` and a
  multiple of `thin`.

- thin:

  thining period for the MCMC sampler. Must be `>= 2`.

- burn:

  burn-in period for the MCMC sampler. Must be in the range
  `1 <= burn < n` and a multiple of `thin`.

- regression:

  logical. Whether to use regression to identify over-dispersion. See
  [`BASiCS_MCMC`](https://rdrr.io/pkg/BASiCS/man/BASiCS_MCMC.html) for
  details.

- params:

  BASiCSParams object to store estimated values in.

- verbose:

  logical. Whether to print progress messages.

- progress:

  logical. Whether to print additional BASiCS progress messages.

- ...:

  Optional parameters passed to
  [`BASiCS_MCMC`](https://rdrr.io/pkg/BASiCS/man/BASiCS_MCMC.html).

## Value

BASiCSParams object containing the estimated parameters.

## Details

This function is just a wrapper around
[`BASiCS_MCMC`](https://rdrr.io/pkg/BASiCS/man/BASiCS_MCMC.html) that
takes the output and converts it to a BASiCSParams object. Either a set
of spike-ins or batch information (or both) must be supplied. If only
batch information is provided there must be at least two batches. See
[`BASiCS_MCMC`](https://rdrr.io/pkg/BASiCS/man/BASiCS_MCMC.html) for
details.

## Examples

``` r
# \donttest{
# Load example data
library(scuttle)
set.seed(1)
sce <- mockSCE()

spike.info <- data.frame(
    Name = rownames(sce)[1:10],
    Input = rnorm(10, 500, 200),
    stringsAsFactors = FALSE
)
params <- BASiCSEstimate(sce[1:100, 1:30], spike.info)
#> Warning: 'BASiCS::newBASiCS_Data' is deprecated.
#> See help("Deprecated")
#> altExp 'spike-ins' is assumed to contain spike-in genes.
#> see help(altExp) for details. 
#> Running with spikes BASiCS sampler (regression case) ... 
#> -----------------------------------------------------
#> MCMC sampler has been started: 20000 iterations to go.
#> -----------------------------------------------------
#> -----------------------------------------------------
#> Iteration 0 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 612.393
#> delta (gene 1): 0.743416
#> phi (cell 1): 0.757726
#> s (cell 1): 1.44699
#> nu (cell 1): 0.416089
#> theta (batch 1): 0.41371
#> betas:    0.8927  -0.0532   0.8013  -1.9110   0.9116   0.9607   0.8638  -1.2077  -0.0110   0.3911  -0.1379  -1.0168
#> 
#> sigma: 0.403606
#> lambda (gene 1): 0.332704
#> epsilon (gene 1): -0.583321
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -4
#> LSdelta (gene 1): -2
#> LSphi: 11
#> LSnu (cell 1): -8.11822
#> LStheta (batch 1): -4
#> -----------------------------------------------------
#> Iteration 20 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 708.912
#> delta (gene 1): 0.627731
#> phi (cell 1): 0.772274
#> s (cell 1): 1.54622
#> nu (cell 1): 0.427917
#> theta (batch 1): 0.309981
#> betas:    0.4507   0.0161   0.1295   0.4678   0.6508   0.6737   0.3515   0.1768  -0.0793  -0.4097  -0.3309  -0.4382
#> 
#> sigma: 0.15265
#> lambda (gene 1): 0.732282
#> epsilon (gene 1): -0.213893
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -4
#> LSdelta (gene 1): -2
#> LSphi: 11
#> LSnu (cell 1): -8.11822
#> LStheta (batch 1): -4
#> -----------------------------------------------------
#> Iteration 40 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 713.338
#> delta (gene 1): 0.724656
#> phi (cell 1): 0.772583
#> s (cell 1): 0.620469
#> nu (cell 1): 0.418365
#> theta (batch 1): 0.366641
#> betas:    1.0819  -0.1418  -0.0693   1.3638   0.6831   0.3289   0.6027  -0.0478   0.4212  -0.5320   0.0834  -0.5715
#> 
#> sigma: 0.152528
#> lambda (gene 1): 0.906016
#> epsilon (gene 1): -0.0976622
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -4
#> LSdelta (gene 1): -2
#> LSphi: 11
#> LSnu (cell 1): -8.11822
#> LStheta (batch 1): -4
#> -----------------------------------------------------
#> Iteration 60 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 637.607
#> delta (gene 1): 0.726841
#> phi (cell 1): 0.77136
#> s (cell 1): 0.544911
#> nu (cell 1): 0.423369
#> theta (batch 1): 0.152667
#> betas:    0.5781  -0.0559   0.6705   1.0208   1.0147   0.0857   0.8605  -0.3467   0.6287  -0.1286  -0.7918   0.2101
#> 
#> sigma: 0.20249
#> lambda (gene 1): 0.983934
#> epsilon (gene 1): 0.0406453
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.9
#> LSdelta (gene 1): -1.9
#> LSphi: 10.9
#> LSnu (cell 1): -8.01822
#> LStheta (batch 1): -3.9
#> -----------------------------------------------------
#> Iteration 80 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 685.889
#> delta (gene 1): 0.618394
#> phi (cell 1): 0.762509
#> s (cell 1): 0.522042
#> nu (cell 1): 0.411238
#> theta (batch 1): 0.149572
#> betas:    0.4904  -0.1300   1.0805   0.7494   0.9214   0.3595   0.5359   0.4569   0.7136  -0.4543  -0.0570   0.1283
#> 
#> sigma: 0.169537
#> lambda (gene 1): 1.93737
#> epsilon (gene 1): -0.0072078
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.9
#> LSdelta (gene 1): -1.9
#> LSphi: 10.9
#> LSnu (cell 1): -8.01822
#> LStheta (batch 1): -3.9
#> -----------------------------------------------------
#> Iteration 100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 563.287
#> delta (gene 1): 0.718131
#> phi (cell 1): 0.760514
#> s (cell 1): 0.407932
#> nu (cell 1): 0.424516
#> theta (batch 1): 0.225828
#> betas:    0.7956  -0.0882   0.6865   0.8374   0.6594   0.8342   0.2556  -0.1392   0.1230  -0.2216  -0.2294  -0.2493
#> 
#> sigma: 0.194888
#> lambda (gene 1): 1.41368
#> epsilon (gene 1): -0.21561
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.8
#> LSdelta (gene 1): -1.8
#> LSphi: 10.8
#> LSnu (cell 1): -7.91822
#> LStheta (batch 1): -3.8
#> -----------------------------------------------------
#> Iteration 120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 797.853
#> delta (gene 1): 0.890289
#> phi (cell 1): 0.782918
#> s (cell 1): 0.972972
#> nu (cell 1): 0.446053
#> theta (batch 1): 0.478069
#> betas:    0.5193  -0.1260   0.0608   1.0022   1.3402  -0.0395   1.2248   0.3314  -0.1100   0.0964  -0.2424   0.2598
#> 
#> sigma: 0.179821
#> lambda (gene 1): 0.603137
#> epsilon (gene 1): 0.0906871
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.8
#> LSdelta (gene 1): -1.8
#> LSphi: 10.8
#> LSnu (cell 1): -7.91822
#> LStheta (batch 1): -3.8
#> -----------------------------------------------------
#> Iteration 140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 735.746
#> delta (gene 1): 0.980585
#> phi (cell 1): 0.801773
#> s (cell 1): 0.404724
#> nu (cell 1): 0.418923
#> theta (batch 1): 0.689914
#> betas:    0.8688  -0.0349   0.4958   0.7103   0.7183   0.6826   0.1785  -0.0599  -0.2749  -0.0929  -0.2423  -0.5017
#> 
#> sigma: 0.186253
#> lambda (gene 1): 0.561452
#> epsilon (gene 1): -0.0976604
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.8
#> LSdelta (gene 1): -1.8
#> LSphi: 10.8
#> LSnu (cell 1): -7.91822
#> LStheta (batch 1): -3.8
#> -----------------------------------------------------
#> Iteration 160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 849.564
#> delta (gene 1): 0.469552
#> phi (cell 1): 0.807331
#> s (cell 1): 2.35875
#> nu (cell 1): 0.416869
#> theta (batch 1): 0.735792
#> betas:    0.6103  -0.0757   0.7574   0.8956   0.5797   1.0035   0.2244   0.0442   0.2758  -0.2695  -0.0394  -0.5081
#> 
#> sigma: 0.23714
#> lambda (gene 1): 0.359041
#> epsilon (gene 1): -0.616668
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.7
#> LSdelta (gene 1): -1.7
#> LSphi: 10.7
#> LSnu (cell 1): -7.81822
#> LStheta (batch 1): -3.7
#> -----------------------------------------------------
#> Iteration 180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 815.994
#> delta (gene 1): 0.810332
#> phi (cell 1): 0.811244
#> s (cell 1): 1.28391
#> nu (cell 1): 0.409256
#> theta (batch 1): 0.340313
#> betas:    0.8094  -0.1115   0.7416   0.8280   0.2467   1.2074  -0.1776   0.1818   0.7682  -0.6071  -0.2893   0.0078
#> 
#> sigma: 0.130827
#> lambda (gene 1): 0.547525
#> epsilon (gene 1): 0.00715975
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.7
#> LSdelta (gene 1): -1.7
#> LSphi: 10.7
#> LSnu (cell 1): -7.81822
#> LStheta (batch 1): -3.7
#> -----------------------------------------------------
#> Iteration 200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 575.594
#> delta (gene 1): 1.36356
#> phi (cell 1): 0.809938
#> s (cell 1): 0.476737
#> nu (cell 1): 0.432895
#> theta (batch 1): 0.13821
#> betas:    0.3955   0.0764   0.7008   0.5727   0.8916   0.3778   0.9189   0.1095  -0.4930   0.2703  -0.8258  -0.9452
#> 
#> sigma: 0.209723
#> lambda (gene 1): 1.9519
#> epsilon (gene 1): 0.486004
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.6
#> LSdelta (gene 1): -1.8
#> LSphi: 10.6
#> LSnu (cell 1): -7.71822
#> LStheta (batch 1): -3.6
#> -----------------------------------------------------
#> Iteration 220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 659.728
#> delta (gene 1): 0.604814
#> phi (cell 1): 0.81435
#> s (cell 1): 0.486223
#> nu (cell 1): 0.413457
#> theta (batch 1): 0.30231
#> betas:    0.1300   0.0509   0.7567   0.5693   1.1919   0.9970   0.1974   0.5449   0.0692  -0.3362  -0.1581  -0.6566
#> 
#> sigma: 0.169631
#> lambda (gene 1): 0.773093
#> epsilon (gene 1): -0.357493
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.6
#> LSdelta (gene 1): -1.8
#> LSphi: 10.6
#> LSnu (cell 1): -7.71822
#> LStheta (batch 1): -3.6
#> -----------------------------------------------------
#> Iteration 240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 732.354
#> delta (gene 1): 1.27007
#> phi (cell 1): 0.832492
#> s (cell 1): 0.553175
#> nu (cell 1): 0.41218
#> theta (batch 1): 0.34972
#> betas:    0.8022  -0.1111   0.2219   0.8808   1.2095   0.7319   0.2099   0.4930   0.1415  -0.5444   0.3313  -0.6329
#> 
#> sigma: 0.16135
#> lambda (gene 1): 1.30039
#> epsilon (gene 1): 0.502133
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.6
#> LSdelta (gene 1): -1.8
#> LSphi: 10.6
#> LSnu (cell 1): -7.71822
#> LStheta (batch 1): -3.6
#> -----------------------------------------------------
#> Iteration 260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 647.068
#> delta (gene 1): 0.807461
#> phi (cell 1): 0.829129
#> s (cell 1): 0.460451
#> nu (cell 1): 0.430204
#> theta (batch 1): 0.253096
#> betas:    1.5201  -0.2056   0.0251   0.9946   1.1511   0.4459   0.2233   0.3148   0.0647  -0.2690   0.1499  -0.7133
#> 
#> sigma: 0.155382
#> lambda (gene 1): 0.569261
#> epsilon (gene 1): -0.0113345
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.5
#> LSdelta (gene 1): -1.7
#> LSphi: 10.5
#> LSnu (cell 1): -7.61822
#> LStheta (batch 1): -3.5
#> -----------------------------------------------------
#> Iteration 280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 639.978
#> delta (gene 1): 0.666937
#> phi (cell 1): 0.853342
#> s (cell 1): 0.553094
#> nu (cell 1): 0.416717
#> theta (batch 1): 0.225468
#> betas:    1.0567  -0.1595   0.2920   0.8114   0.7284   0.2207   0.6833   0.4518  -0.0570  -0.2470  -0.0401  -0.0734
#> 
#> sigma: 0.138869
#> lambda (gene 1): 1.10543
#> epsilon (gene 1): -0.147818
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.5
#> LSdelta (gene 1): -1.7
#> LSphi: 10.5
#> LSnu (cell 1): -7.61822
#> LStheta (batch 1): -3.5
#> -----------------------------------------------------
#> Iteration 300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 675.92
#> delta (gene 1): 0.604439
#> phi (cell 1): 0.837323
#> s (cell 1): 1.3003
#> nu (cell 1): 0.426241
#> theta (batch 1): 0.24225
#> betas:    0.9732  -0.0926   0.2612  -0.1928   0.4569   1.1516   0.2826   0.0395  -0.0136   0.1290  -0.6122  -0.2708
#> 
#> sigma: 0.201467
#> lambda (gene 1): 0.63684
#> epsilon (gene 1): -0.405739
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.4
#> LSdelta (gene 1): -1.6
#> LSphi: 10.4
#> LSnu (cell 1): -7.51822
#> LStheta (batch 1): -3.4
#> -----------------------------------------------------
#> Iteration 320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 655.698
#> delta (gene 1): 0.876922
#> phi (cell 1): 0.82786
#> s (cell 1): 0.31607
#> nu (cell 1): 0.415104
#> theta (batch 1): 0.350383
#> betas:    0.3184  -0.0387   0.0643   1.1955   0.4183   0.8956   0.6721   0.0579   0.4675  -0.0547  -0.4268  -0.0155
#> 
#> sigma: 0.124401
#> lambda (gene 1): 0.380198
#> epsilon (gene 1): -0.0825703
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.4
#> LSdelta (gene 1): -1.6
#> LSphi: 10.4
#> LSnu (cell 1): -7.51822
#> LStheta (batch 1): -3.4
#> -----------------------------------------------------
#> Iteration 340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 867.198
#> delta (gene 1): 1.03997
#> phi (cell 1): 0.818188
#> s (cell 1): 0.32805
#> nu (cell 1): 0.42971
#> theta (batch 1): 0.146339
#> betas:    1.0133  -0.1006   0.0149   0.6791   0.6707   0.8055   0.1387   0.4708  -0.1521  -0.1914  -0.0678  -0.5161
#> 
#> sigma: 0.132267
#> lambda (gene 1): 0.752804
#> epsilon (gene 1): 0.108233
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.4
#> LSdelta (gene 1): -1.6
#> LSphi: 10.4
#> LSnu (cell 1): -7.51822
#> LStheta (batch 1): -3.4
#> -----------------------------------------------------
#> Iteration 360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 713.789
#> delta (gene 1): 0.472612
#> phi (cell 1): 0.832207
#> s (cell 1): 0.517525
#> nu (cell 1): 0.421938
#> theta (batch 1): 0.262531
#> betas:    0.1679   0.0257   0.1376   1.1785   0.1898   1.3269   0.3361   0.2596   0.4829  -0.3532  -0.3635  -0.4714
#> 
#> sigma: 0.176775
#> lambda (gene 1): 0.920756
#> epsilon (gene 1): -0.701453
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.3
#> LSdelta (gene 1): -1.5
#> LSphi: 10.3
#> LSnu (cell 1): -7.41822
#> LStheta (batch 1): -3.3
#> -----------------------------------------------------
#> Iteration 380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 736.781
#> delta (gene 1): 0.666058
#> phi (cell 1): 0.821979
#> s (cell 1): 0.372263
#> nu (cell 1): 0.424378
#> theta (batch 1): 0.114052
#> betas:    1.3384  -0.0973  -0.2927   0.1801   1.0855   0.1801   0.4275  -0.0618  -0.0486  -0.1928  -0.5430  -0.4732
#> 
#> sigma: 0.190045
#> lambda (gene 1): 0.315672
#> epsilon (gene 1): -0.356833
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.3
#> LSdelta (gene 1): -1.5
#> LSphi: 10.3
#> LSnu (cell 1): -7.41822
#> LStheta (batch 1): -3.3
#> -----------------------------------------------------
#> Iteration 400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 760.011
#> delta (gene 1): 0.535978
#> phi (cell 1): 0.801485
#> s (cell 1): 0.30151
#> nu (cell 1): 0.431843
#> theta (batch 1): 0.201115
#> betas:    0.7275  -0.0449   0.4524   0.3625   0.8253   0.8039   0.2767   0.1012  -0.0473   0.0478  -0.5399  -0.4009
#> 
#> sigma: 0.165869
#> lambda (gene 1): 1.10833
#> epsilon (gene 1): -0.55017
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.2
#> LSdelta (gene 1): -1.4
#> LSphi: 10.2
#> LSnu (cell 1): -7.31822
#> LStheta (batch 1): -3.2
#> -----------------------------------------------------
#> Iteration 420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 797.011
#> delta (gene 1): 0.586877
#> phi (cell 1): 0.800076
#> s (cell 1): 0.283104
#> nu (cell 1): 0.403224
#> theta (batch 1): 0.251445
#> betas:    1.4702  -0.1754  -0.1947   0.3712   0.9484   0.3178   0.1274  -0.0425   0.6061  -0.6747  -0.2841  -0.1055
#> 
#> sigma: 0.190994
#> lambda (gene 1): 0.524659
#> epsilon (gene 1): -0.296731
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.2
#> LSdelta (gene 1): -1.4
#> LSphi: 10.2
#> LSnu (cell 1): -7.31822
#> LStheta (batch 1): -3.2
#> -----------------------------------------------------
#> Iteration 440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1032.46
#> delta (gene 1): 0.639972
#> phi (cell 1): 0.764849
#> s (cell 1): 1.14933
#> nu (cell 1): 0.426301
#> theta (batch 1): 0.397822
#> betas:    0.4003  -0.0649   0.7473   0.5183   0.5497   0.9688   0.8513   0.0520   0.4991  -0.2633  -0.2284  -0.0539
#> 
#> sigma: 0.178163
#> lambda (gene 1): 0.607369
#> epsilon (gene 1): -0.308058
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.2
#> LSdelta (gene 1): -1.4
#> LSphi: 10.2
#> LSnu (cell 1): -7.31822
#> LStheta (batch 1): -3.2
#> -----------------------------------------------------
#> Iteration 460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 584.2
#> delta (gene 1): 0.646068
#> phi (cell 1): 0.806378
#> s (cell 1): 0.932742
#> nu (cell 1): 0.415482
#> theta (batch 1): 0.488112
#> betas:    0.6290  -0.0993   0.1958  -0.2268   1.1604   0.7165   0.3500   0.6171   0.1234   0.1905  -0.2630  -0.2491
#> 
#> sigma: 0.145881
#> lambda (gene 1): 1.66148
#> epsilon (gene 1): -0.280509
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.1
#> LSdelta (gene 1): -1.3
#> LSphi: 10.1
#> LSnu (cell 1): -7.21822
#> LStheta (batch 1): -3.1
#> -----------------------------------------------------
#> Iteration 480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 795.977
#> delta (gene 1): 0.638956
#> phi (cell 1): 0.81262
#> s (cell 1): 0.2123
#> nu (cell 1): 0.418381
#> theta (batch 1): 0.228252
#> betas:    0.8413  -0.0736   0.0904   0.3983   0.6103   0.6890   1.0004   0.0987   0.1261  -0.0640  -0.3314  -0.4301
#> 
#> sigma: 0.17319
#> lambda (gene 1): 0.884156
#> epsilon (gene 1): -0.25634
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3.1
#> LSdelta (gene 1): -1.3
#> LSphi: 10.1
#> LSnu (cell 1): -7.21822
#> LStheta (batch 1): -3.1
#> -----------------------------------------------------
#> Iteration 500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 648.616
#> delta (gene 1): 0.776702
#> phi (cell 1): 0.81141
#> s (cell 1): 0.58267
#> nu (cell 1): 0.4202
#> theta (batch 1): 0.174408
#> betas:    0.5528  -0.0604   0.4497   0.7047   0.5900   0.9806   0.3068   0.5095   0.1271  -0.3055   0.0225  -0.5205
#> 
#> sigma: 0.239332
#> lambda (gene 1): 0.907951
#> epsilon (gene 1): -0.146181
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3
#> LSdelta (gene 1): -1.2
#> LSphi: 10
#> LSnu (cell 1): -7.11822
#> LStheta (batch 1): -3
#> -----------------------------------------------------
#> Iteration 520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 827.89
#> delta (gene 1): 0.632356
#> phi (cell 1): 0.801272
#> s (cell 1): 0.53867
#> nu (cell 1): 0.424486
#> theta (batch 1): 0.283708
#> betas:    0.6219   0.0057   0.2170   0.3281   1.2539   0.6854   0.5134   0.2618  -0.2961  -0.1285  -0.3374  -0.4431
#> 
#> sigma: 0.188041
#> lambda (gene 1): 1.0057
#> epsilon (gene 1): -0.462772
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3
#> LSdelta (gene 1): -1.2
#> LSphi: 10
#> LSnu (cell 1): -7.11822
#> LStheta (batch 1): -3
#> -----------------------------------------------------
#> Iteration 540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 880.439
#> delta (gene 1): 1.51701
#> phi (cell 1): 0.792955
#> s (cell 1): 1.05777
#> nu (cell 1): 0.411961
#> theta (batch 1): 0.197289
#> betas:    0.2561  -0.0943   0.0698   1.0777   0.2390   1.1243   0.8561   0.5585   0.1516   0.0286   0.0278  -0.0103
#> 
#> sigma: 0.119578
#> lambda (gene 1): 0.507306
#> epsilon (gene 1): 0.500132
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -3
#> LSdelta (gene 1): -1.2
#> LSphi: 10
#> LSnu (cell 1): -7.11822
#> LStheta (batch 1): -3
#> -----------------------------------------------------
#> Iteration 560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 553.995
#> delta (gene 1): 1.02851
#> phi (cell 1): 0.784845
#> s (cell 1): 0.405615
#> nu (cell 1): 0.408771
#> theta (batch 1): 0.327566
#> betas:    0.4563  -0.0400  -0.2376   1.2150   0.7709   0.6890   0.0596   0.5470   0.3958  -0.7887   0.4626  -0.7534
#> 
#> sigma: 0.186469
#> lambda (gene 1): 2.1682
#> epsilon (gene 1): 0.153687
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.9
#> LSdelta (gene 1): -1.3
#> LSphi: 9.9
#> LSnu (cell 1): -7.01822
#> LStheta (batch 1): -2.9
#> -----------------------------------------------------
#> Iteration 580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 632.676
#> delta (gene 1): 0.508431
#> phi (cell 1): 0.810843
#> s (cell 1): 0.178114
#> nu (cell 1): 0.408561
#> theta (batch 1): 0.340004
#> betas:    0.6527  -0.0053   0.3595   0.3598   1.0824   0.6626   0.4551  -0.0985   0.2627  -0.7109  -0.0458  -0.7096
#> 
#> sigma: 0.139997
#> lambda (gene 1): 0.4616
#> epsilon (gene 1): -0.470895
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.9
#> LSdelta (gene 1): -1.3
#> LSphi: 9.9
#> LSnu (cell 1): -7.01822
#> LStheta (batch 1): -2.9
#> -----------------------------------------------------
#> Iteration 600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 931.793
#> delta (gene 1): 0.598007
#> phi (cell 1): 0.79078
#> s (cell 1): 0.460381
#> nu (cell 1): 0.420494
#> theta (batch 1): 0.294981
#> betas:    0.8558   0.0060   0.6834   1.6290  -0.2696   1.0643   0.4020  -0.2071  -0.1039  -0.1532  -0.9213  -0.4225
#> 
#> sigma: 0.191649
#> lambda (gene 1): 1.30539
#> epsilon (gene 1): -0.327858
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.8
#> LSdelta (gene 1): -1.2
#> LSphi: 9.8
#> LSnu (cell 1): -6.91822
#> LStheta (batch 1): -2.8
#> -----------------------------------------------------
#> Iteration 620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 663.973
#> delta (gene 1): 0.886481
#> phi (cell 1): 0.797743
#> s (cell 1): 0.366962
#> nu (cell 1): 0.412171
#> theta (batch 1): 0.237795
#> betas:    0.7508  -0.0271  -0.4204   1.1378   0.8277   0.2386   0.2935   0.4584   0.0124  -0.7725   0.0609  -0.4495
#> 
#> sigma: 0.156683
#> lambda (gene 1): 1.02235
#> epsilon (gene 1): 0.04849
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.8
#> LSdelta (gene 1): -1.2
#> LSphi: 9.8
#> LSnu (cell 1): -6.91822
#> LStheta (batch 1): -2.8
#> -----------------------------------------------------
#> Iteration 640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 684.189
#> delta (gene 1): 0.851437
#> phi (cell 1): 0.746077
#> s (cell 1): 0.471177
#> nu (cell 1): 0.415212
#> theta (batch 1): 0.134536
#> betas:    1.1088  -0.0175   0.7671   0.6808   0.3287   0.6300   0.4514   0.0478  -0.2968  -0.6187  -0.3521  -0.6459
#> 
#> sigma: 0.150178
#> lambda (gene 1): 1.08392
#> epsilon (gene 1): 0.0874749
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.8
#> LSdelta (gene 1): -1.2
#> LSphi: 9.8
#> LSnu (cell 1): -6.91822
#> LStheta (batch 1): -2.8
#> -----------------------------------------------------
#> Iteration 660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 619.803
#> delta (gene 1): 0.619583
#> phi (cell 1): 0.747034
#> s (cell 1): 0.516548
#> nu (cell 1): 0.40947
#> theta (batch 1): 0.0528125
#> betas:    0.3487  -0.0092   0.6150   0.7591   0.9397   1.1313  -0.2640   0.4216   0.2442  -0.4164  -0.5328  -0.1109
#> 
#> sigma: 0.184606
#> lambda (gene 1): 0.496853
#> epsilon (gene 1): -0.235391
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.7
#> LSdelta (gene 1): -1.3
#> LSphi: 9.7
#> LSnu (cell 1): -6.81822
#> LStheta (batch 1): -2.7
#> -----------------------------------------------------
#> Iteration 680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 719.477
#> delta (gene 1): 0.816503
#> phi (cell 1): 0.780245
#> s (cell 1): 0.327766
#> nu (cell 1): 0.420606
#> theta (batch 1): 0.115936
#> betas:    1.2695  -0.1285   0.8128   0.4609   0.3874   0.7961   0.2682  -0.4213   0.1891  -0.1427  -0.6513   0.0576
#> 
#> sigma: 0.156466
#> lambda (gene 1): 0.568002
#> epsilon (gene 1): -0.107946
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.7
#> LSdelta (gene 1): -1.3
#> LSphi: 9.7
#> LSnu (cell 1): -6.81822
#> LStheta (batch 1): -2.7
#> -----------------------------------------------------
#> Iteration 700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 604.417
#> delta (gene 1): 0.717157
#> phi (cell 1): 0.761547
#> s (cell 1): 0.337213
#> nu (cell 1): 0.426637
#> theta (batch 1): 0.187597
#> betas:    0.8709   0.0536   0.4402   0.5335   0.5962   0.2676   0.1299   0.2932  -0.5380  -0.7042  -0.1594  -1.0690
#> 
#> sigma: 0.136749
#> lambda (gene 1): 0.640002
#> epsilon (gene 1): -0.136283
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.6
#> LSdelta (gene 1): -1.4
#> LSphi: 9.6
#> LSnu (cell 1): -6.71822
#> LStheta (batch 1): -2.6
#> -----------------------------------------------------
#> Iteration 720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 613.447
#> delta (gene 1): 0.402623
#> phi (cell 1): 0.710469
#> s (cell 1): 0.677531
#> nu (cell 1): 0.439273
#> theta (batch 1): 0.218764
#> betas:    0.7581  -0.0063   0.5707   0.8501   0.8668   0.3730   0.1744   0.3671  -0.1347  -0.4393  -0.0967  -0.9864
#> 
#> sigma: 0.175653
#> lambda (gene 1): 0.617058
#> epsilon (gene 1): -0.830177
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.6
#> LSdelta (gene 1): -1.4
#> LSphi: 9.6
#> LSnu (cell 1): -6.71822
#> LStheta (batch 1): -2.6
#> -----------------------------------------------------
#> Iteration 740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 703.825
#> delta (gene 1): 0.600101
#> phi (cell 1): 0.730204
#> s (cell 1): 0.496036
#> nu (cell 1): 0.429001
#> theta (batch 1): 0.095981
#> betas:    0.7504  -0.0573   0.7018   0.8750   0.4853   0.8333   0.1889   0.1142   0.0929  -0.1654  -0.6862   0.0147
#> 
#> sigma: 0.123243
#> lambda (gene 1): 0.363062
#> epsilon (gene 1): -0.303909
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.6
#> LSdelta (gene 1): -1.4
#> LSphi: 9.6
#> LSnu (cell 1): -6.71822
#> LStheta (batch 1): -2.6
#> -----------------------------------------------------
#> Iteration 760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 832.839
#> delta (gene 1): 0.689339
#> phi (cell 1): 0.722997
#> s (cell 1): 0.35456
#> nu (cell 1): 0.419765
#> theta (batch 1): 0.0347035
#> betas:    0.4402  -0.0168   0.2834   0.4494   0.9795   0.8666   0.1643   0.4259  -0.2491   0.0835  -0.3818  -0.6677
#> 
#> sigma: 0.154
#> lambda (gene 1): 0.764227
#> epsilon (gene 1): -0.0977833
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.3
#> LSphi: 9.5
#> LSnu (cell 1): -6.61822
#> LStheta (batch 1): -2.5
#> -----------------------------------------------------
#> Iteration 780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 807.264
#> delta (gene 1): 0.75734
#> phi (cell 1): 0.688066
#> s (cell 1): 0.57054
#> nu (cell 1): 0.425891
#> theta (batch 1): 0.127587
#> betas:    0.8475   0.0022   0.0081   0.4986   0.8734   0.2697   0.4200   0.2678  -0.3586  -0.2617  -0.3500  -0.8854
#> 
#> sigma: 0.199738
#> lambda (gene 1): 0.717843
#> epsilon (gene 1): -0.0791421
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.3
#> LSphi: 9.5
#> LSnu (cell 1): -6.61822
#> LStheta (batch 1): -2.5
#> -----------------------------------------------------
#> Iteration 800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 788.281
#> delta (gene 1): 0.489018
#> phi (cell 1): 0.738075
#> s (cell 1): 0.398237
#> nu (cell 1): 0.42839
#> theta (batch 1): 0.129938
#> betas:    1.0904  -0.0596   0.5966  -0.1029   0.8478   0.4256   0.6766  -0.2720  -0.0559   0.0724  -0.8450  -0.2829
#> 
#> sigma: 0.151551
#> lambda (gene 1): 0.846606
#> epsilon (gene 1): -0.554855
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 9.4
#> LSnu (cell 1): -6.51822
#> LStheta (batch 1): -2.4
#> -----------------------------------------------------
#> Iteration 820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 738.774
#> delta (gene 1): 0.58634
#> phi (cell 1): 0.67641
#> s (cell 1): 0.465747
#> nu (cell 1): 0.424523
#> theta (batch 1): 0.179922
#> betas:    1.1874  -0.1319   0.6333   0.4089   0.7503   0.8865   0.0589   0.4373  -0.2225  -0.1759  -0.4422   0.0644
#> 
#> sigma: 0.127953
#> lambda (gene 1): 3.22048
#> epsilon (gene 1): -0.258528
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 9.4
#> LSnu (cell 1): -6.51822
#> LStheta (batch 1): -2.4
#> -----------------------------------------------------
#> Iteration 840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 830.355
#> delta (gene 1): 0.528007
#> phi (cell 1): 0.687979
#> s (cell 1): 0.496224
#> nu (cell 1): 0.430086
#> theta (batch 1): 0.298198
#> betas:    0.2802   0.0281   0.0636   1.0447   1.1983   0.3193   1.0277  -0.1399   0.0194   0.2352  -0.9752  -0.1554
#> 
#> sigma: 0.226076
#> lambda (gene 1): 0.885336
#> epsilon (gene 1): -0.40377
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 9.4
#> LSnu (cell 1): -6.51822
#> LStheta (batch 1): -2.4
#> -----------------------------------------------------
#> Iteration 860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 685.662
#> delta (gene 1): 0.731814
#> phi (cell 1): 0.648184
#> s (cell 1): 0.331016
#> nu (cell 1): 0.427271
#> theta (batch 1): 0.156761
#> betas:    0.4096  -0.0399   0.1325   1.0910   0.5301   1.2395   0.2842   0.1427  -0.0615   0.0272  -0.3229  -0.0779
#> 
#> sigma: 0.179701
#> lambda (gene 1): 1.15507
#> epsilon (gene 1): -0.27958
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.1
#> LSphi: 9.3
#> LSnu (cell 1): -6.41822
#> LStheta (batch 1): -2.3
#> -----------------------------------------------------
#> Iteration 880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 673.96
#> delta (gene 1): 0.678326
#> phi (cell 1): 0.641687
#> s (cell 1): 0.417818
#> nu (cell 1): 0.433908
#> theta (batch 1): 0.236148
#> betas:    0.9943  -0.1251   0.5496   0.3371   1.1194   0.3360   0.5849  -0.1528   0.3182  -0.2643  -0.3462  -0.1656
#> 
#> sigma: 0.141664
#> lambda (gene 1): 0.873194
#> epsilon (gene 1): -0.278083
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.1
#> LSphi: 9.3
#> LSnu (cell 1): -6.41822
#> LStheta (batch 1): -2.3
#> -----------------------------------------------------
#> Iteration 900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 787.61
#> delta (gene 1): 0.521418
#> phi (cell 1): 0.714554
#> s (cell 1): 0.255872
#> nu (cell 1): 0.424335
#> theta (batch 1): 0.213541
#> betas:    0.7380  -0.1173   0.7957   0.4089   1.0680   0.5434   0.7736   0.3666  -0.1830  -0.1579  -0.0637  -0.0135
#> 
#> sigma: 0.16903
#> lambda (gene 1): 0.525802
#> epsilon (gene 1): -0.369041
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1
#> LSphi: 9.2
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -2.2
#> -----------------------------------------------------
#> Iteration 920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 537.921
#> delta (gene 1): 0.537307
#> phi (cell 1): 0.716713
#> s (cell 1): 0.289781
#> nu (cell 1): 0.427002
#> theta (batch 1): 0.17298
#> betas:    1.1756  -0.1513   0.6975   0.1739   0.5249   1.3194  -0.2577   0.6120   0.0370  -0.1496  -0.7317   0.3980
#> 
#> sigma: 0.232639
#> lambda (gene 1): 0.654818
#> epsilon (gene 1): -0.275992
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1
#> LSphi: 9.2
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -2.2
#> -----------------------------------------------------
#> Iteration 940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1079.57
#> delta (gene 1): 0.967165
#> phi (cell 1): 0.756475
#> s (cell 1): 0.386714
#> nu (cell 1): 0.402175
#> theta (batch 1): 0.300527
#> betas:    0.8387  -0.1039   0.6844   0.5276   0.6850   0.7592   0.2546   0.4344   0.0357  -0.0821  -0.3901  -0.0633
#> 
#> sigma: 0.143005
#> lambda (gene 1): 1.26685
#> epsilon (gene 1): 0.168212
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1
#> LSphi: 9.2
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -2.2
#> -----------------------------------------------------
#> Iteration 960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 779.952
#> delta (gene 1): 0.697679
#> phi (cell 1): 0.705303
#> s (cell 1): 0.805381
#> nu (cell 1): 0.438588
#> theta (batch 1): 0.445646
#> betas:    0.8170  -0.1303   0.6221   0.3136   0.7740   0.9440   0.0692   1.0060  -0.1915  -0.3059   0.0737  -0.3613
#> 
#> sigma: 0.152776
#> lambda (gene 1): 1.19028
#> epsilon (gene 1): -0.167072
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.1
#> LSphi: 9.1
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -2.1
#> -----------------------------------------------------
#> Iteration 980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 610.574
#> delta (gene 1): 0.653536
#> phi (cell 1): 0.741158
#> s (cell 1): 0.197299
#> nu (cell 1): 0.432135
#> theta (batch 1): 0.564678
#> betas:    0.2708   0.0730   0.7979   0.8121   0.7544   0.5069   0.4200   0.3131  -0.2044  -0.4321  -0.5263  -0.6141
#> 
#> sigma: 0.132476
#> lambda (gene 1): 0.488464
#> epsilon (gene 1): -0.216342
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.1
#> LSphi: 9.1
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -2.1
#> -----------------------------------------------------
#> Iteration 1000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 663.244
#> delta (gene 1): 0.545635
#> phi (cell 1): 0.68145
#> s (cell 1): 1.69347
#> nu (cell 1): 0.428388
#> theta (batch 1): 0.477622
#> betas:    1.2419  -0.1729   0.6332   1.1355   0.5225   0.4724   0.2580   0.7642  -0.7842   0.2009  -0.2500  -0.1451
#> 
#> sigma: 0.212931
#> lambda (gene 1): 0.725376
#> epsilon (gene 1): -0.368814
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.2
#> LSphi: 9
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -2
#> -----------------------------------------------------
#> Iteration 1020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 702.356
#> delta (gene 1): 0.885672
#> phi (cell 1): 0.656077
#> s (cell 1): 0.450988
#> nu (cell 1): 0.422775
#> theta (batch 1): 0.106981
#> betas:    0.8141  -0.1402   0.5123   0.7630   0.8367   0.4667   0.5941   0.0717   0.1138   0.0608  -0.3756  -0.0661
#> 
#> sigma: 0.116066
#> lambda (gene 1): 1.82948
#> epsilon (gene 1): 0.0447133
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.2
#> LSphi: 9
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -2
#> -----------------------------------------------------
#> Iteration 1040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 914.22
#> delta (gene 1): 0.68757
#> phi (cell 1): 0.651759
#> s (cell 1): 0.452876
#> nu (cell 1): 0.427612
#> theta (batch 1): 0.0145945
#> betas:    0.6662  -0.1092   0.8910   1.2868   0.3977   0.4389   0.5399   0.4241  -0.0730   0.1625  -0.5544   0.3225
#> 
#> sigma: 0.131622
#> lambda (gene 1): 0.928615
#> epsilon (gene 1): -0.286783
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.2
#> LSphi: 9
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -2
#> -----------------------------------------------------
#> Iteration 1060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 678.072
#> delta (gene 1): 0.549119
#> phi (cell 1): 0.663242
#> s (cell 1): 0.409934
#> nu (cell 1): 0.424727
#> theta (batch 1): 0.006108
#> betas:    0.7788   0.0498   0.4761   0.1200   0.6444   0.6841  -0.0054   0.2518  -0.5364  -0.2764  -0.5835  -1.0292
#> 
#> sigma: 0.171529
#> lambda (gene 1): 0.282528
#> epsilon (gene 1): -0.472393
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.1
#> LSphi: 8.9
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.9
#> -----------------------------------------------------
#> Iteration 1080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 874.173
#> delta (gene 1): 0.571287
#> phi (cell 1): 0.626282
#> s (cell 1): 0.381455
#> nu (cell 1): 0.417287
#> theta (batch 1): 0.0124832
#> betas:    1.1566  -0.2190   0.7742   0.5578   0.6670   0.8444   0.1452   0.5202   0.1369  -0.2014   0.0345   0.4444
#> 
#> sigma: 0.167778
#> lambda (gene 1): 2.01693
#> epsilon (gene 1): -0.437821
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.1
#> LSphi: 8.9
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.9
#> -----------------------------------------------------
#> Iteration 1100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 767.156
#> delta (gene 1): 0.733467
#> phi (cell 1): 0.632049
#> s (cell 1): 0.858189
#> nu (cell 1): 0.425863
#> theta (batch 1): 0.126462
#> betas:    0.9483  -0.0449   0.5822   0.6511   1.3910  -0.2976   0.6578   0.0150  -0.3469   0.4143  -1.1905  -0.2615
#> 
#> sigma: 0.187103
#> lambda (gene 1): 1.30044
#> epsilon (gene 1): -0.00642973
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2
#> LSdelta (gene 1): -1
#> LSphi: 8.8
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.8
#> -----------------------------------------------------
#> Iteration 1120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 621.593
#> delta (gene 1): 0.830716
#> phi (cell 1): 0.646653
#> s (cell 1): 0.423775
#> nu (cell 1): 0.423568
#> theta (batch 1): 0.0834497
#> betas:    0.7945  -0.1595   1.3688   0.9292   0.6535   0.9528   0.3051   0.6463   0.1316   0.1720  -0.6054   0.4037
#> 
#> sigma: 0.142972
#> lambda (gene 1): 1.44714
#> epsilon (gene 1): 0.137464
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2
#> LSdelta (gene 1): -1
#> LSphi: 8.8
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.8
#> -----------------------------------------------------
#> Iteration 1140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 793.168
#> delta (gene 1): 0.854839
#> phi (cell 1): 0.683607
#> s (cell 1): 0.686585
#> nu (cell 1): 0.412168
#> theta (batch 1): 0.15167
#> betas:    0.8392  -0.0684   0.3868   1.2738   0.4617   0.4112   0.8668   0.0491  -0.1187   0.3628  -0.7900  -0.2030
#> 
#> sigma: 0.206803
#> lambda (gene 1): 1.39393
#> epsilon (gene 1): -0.0444602
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2
#> LSdelta (gene 1): -1
#> LSphi: 8.8
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.8
#> -----------------------------------------------------
#> Iteration 1160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 804.472
#> delta (gene 1): 0.697458
#> phi (cell 1): 0.709274
#> s (cell 1): 0.437922
#> nu (cell 1): 0.424615
#> theta (batch 1): 0.054089
#> betas:    0.6975  -0.0923   0.9728   0.5859   0.5449   0.8225   0.2736   0.2855   0.2295  -0.2193  -0.0866  -0.2944
#> 
#> sigma: 0.13888
#> lambda (gene 1): 3.32549
#> epsilon (gene 1): -0.240774
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.1
#> LSphi: 8.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.7
#> -----------------------------------------------------
#> Iteration 1180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 776.084
#> delta (gene 1): 0.629867
#> phi (cell 1): 0.693223
#> s (cell 1): 0.515936
#> nu (cell 1): 0.435242
#> theta (batch 1): 0.562321
#> betas:    1.1133  -0.1419   0.1836   0.3408   0.7949   0.8281  -0.0406   0.3287   0.4249  -0.9232   0.2344  -0.3340
#> 
#> sigma: 0.147938
#> lambda (gene 1): 1.26493
#> epsilon (gene 1): -0.104356
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.1
#> LSphi: 8.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.7
#> -----------------------------------------------------
#> Iteration 1200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 736.946
#> delta (gene 1): 0.570635
#> phi (cell 1): 0.722897
#> s (cell 1): 0.675797
#> nu (cell 1): 0.430976
#> theta (batch 1): 0.451448
#> betas:    1.2230  -0.2017   0.3999   0.8816   0.4857   1.0746  -0.1528   0.5318   0.3859  -0.4172  -0.1372   0.1279
#> 
#> sigma: 0.152779
#> lambda (gene 1): 1.13874
#> epsilon (gene 1): -0.242877
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.2
#> LSphi: 8.6
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.8
#> -----------------------------------------------------
#> Iteration 1220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 754.648
#> delta (gene 1): 0.560472
#> phi (cell 1): 0.720945
#> s (cell 1): 0.560709
#> nu (cell 1): 0.428232
#> theta (batch 1): 0.316404
#> betas:    1.1500  -0.0832   0.5059   0.3107   0.9029   0.5146   0.2197   0.0500   0.1598  -0.0191  -0.8671  -0.1969
#> 
#> sigma: 0.160902
#> lambda (gene 1): 2.70132
#> epsilon (gene 1): -0.372797
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.2
#> LSphi: 8.6
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.8
#> -----------------------------------------------------
#> Iteration 1240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 778.061
#> delta (gene 1): 0.951134
#> phi (cell 1): 0.724203
#> s (cell 1): 0.483046
#> nu (cell 1): 0.424833
#> theta (batch 1): 0.435967
#> betas:    0.6153  -0.0636   0.5644   1.3046   0.6351   1.0517   0.3143   0.3932   0.0860  -0.0545  -0.4018  -0.0307
#> 
#> sigma: 0.169676
#> lambda (gene 1): 0.972218
#> epsilon (gene 1): 0.0895024
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.2
#> LSphi: 8.6
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.8
#> -----------------------------------------------------
#> Iteration 1260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 642.979
#> delta (gene 1): 0.617027
#> phi (cell 1): 0.807762
#> s (cell 1): 3.63308
#> nu (cell 1): 0.418652
#> theta (batch 1): 0.43119
#> betas:    0.8452  -0.0508   0.6978   0.3671   1.0121   0.4820   0.2436   0.0575   0.2898  -0.5366  -0.1860  -0.8013
#> 
#> sigma: 0.151264
#> lambda (gene 1): 1.23261
#> epsilon (gene 1): -0.285548
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 8.5
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.9
#> -----------------------------------------------------
#> Iteration 1280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 603.818
#> delta (gene 1): 0.764581
#> phi (cell 1): 0.804347
#> s (cell 1): 1.05488
#> nu (cell 1): 0.411662
#> theta (batch 1): 0.737403
#> betas:    1.0789  -0.1377   0.4705   0.1983   1.1333   0.4710   0.3336   0.1590   0.0576  -0.0601  -0.3171  -0.0280
#> 
#> sigma: 0.115097
#> lambda (gene 1): 1.00615
#> epsilon (gene 1): -0.217319
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 8.5
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.9
#> -----------------------------------------------------
#> Iteration 1300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 812.007
#> delta (gene 1): 0.565567
#> phi (cell 1): 0.701164
#> s (cell 1): 0.5006
#> nu (cell 1): 0.427138
#> theta (batch 1): 0.373973
#> betas:    0.6830  -0.1047   0.6831   1.4016   0.6222   0.5019  -0.0416   0.7757   0.1215  -0.3311   0.3396  -0.5022
#> 
#> sigma: 0.182196
#> lambda (gene 1): 0.463853
#> epsilon (gene 1): -0.524647
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 8.4
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.8
#> -----------------------------------------------------
#> Iteration 1320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 677.499
#> delta (gene 1): 0.676265
#> phi (cell 1): 0.719676
#> s (cell 1): 0.549148
#> nu (cell 1): 0.429797
#> theta (batch 1): 0.248768
#> betas:    0.8147  -0.0625   0.5522   0.3302   1.0598   0.4410   0.4397   0.2168  -0.3065  -0.0445  -0.5412  -0.3555
#> 
#> sigma: 0.132604
#> lambda (gene 1): 0.775571
#> epsilon (gene 1): -0.0328217
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 8.4
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.8
#> -----------------------------------------------------
#> Iteration 1340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 833.508
#> delta (gene 1): 0.998144
#> phi (cell 1): 0.716289
#> s (cell 1): 0.491407
#> nu (cell 1): 0.405501
#> theta (batch 1): 0.357384
#> betas:    0.8266  -0.2118   0.7581   0.6100   0.6405   0.8152   0.6945   0.3099   0.2025   0.2973  -0.4207   0.8661
#> 
#> sigma: 0.170371
#> lambda (gene 1): 1.84544
#> epsilon (gene 1): 0.160688
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 8.4
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.8
#> -----------------------------------------------------
#> Iteration 1360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1046.44
#> delta (gene 1): 0.604342
#> phi (cell 1): 0.717365
#> s (cell 1): 0.263729
#> nu (cell 1): 0.430588
#> theta (batch 1): 0.320478
#> betas:    1.5887  -0.2518   0.2001   0.5614   0.4293   0.5321   0.4351   0.2112   0.1463   0.0127  -0.6515   0.3242
#> 
#> sigma: 0.151914
#> lambda (gene 1): 2.80238
#> epsilon (gene 1): -0.156087
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.3
#> LSphi: 8.3
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.7
#> -----------------------------------------------------
#> Iteration 1380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 568.716
#> delta (gene 1): 0.799811
#> phi (cell 1): 0.738973
#> s (cell 1): 0.541375
#> nu (cell 1): 0.414271
#> theta (batch 1): 0.10862
#> betas:    0.9737  -0.2365   0.4949   0.8014   1.0662   0.2959   0.9792   0.0392   0.4992  -0.2806   0.3365   0.0971
#> 
#> sigma: 0.196339
#> lambda (gene 1): 0.661097
#> epsilon (gene 1): 0.00936918
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.3
#> LSphi: 8.3
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.7
#> -----------------------------------------------------
#> Iteration 1400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 679.612
#> delta (gene 1): 0.82492
#> phi (cell 1): 0.718328
#> s (cell 1): 0.549443
#> nu (cell 1): 0.426938
#> theta (batch 1): 0.30458
#> betas:    0.2160   0.1167   0.3426   1.1196   0.4321   0.8878   0.6535  -0.2866   0.2543  -0.6064  -0.7158  -0.4168
#> 
#> sigma: 0.137629
#> lambda (gene 1): 2.01123
#> epsilon (gene 1): 0.104146
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 8.2
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.6
#> -----------------------------------------------------
#> Iteration 1420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 652.73
#> delta (gene 1): 0.658058
#> phi (cell 1): 0.707533
#> s (cell 1): 1.59535
#> nu (cell 1): 0.415131
#> theta (batch 1): 0.282955
#> betas:    0.3166  -0.0550   1.2163  -0.2062   1.0306   0.9680   0.5391   0.1930   0.0799  -0.0382  -0.0876  -0.2815
#> 
#> sigma: 0.190694
#> lambda (gene 1): 0.957058
#> epsilon (gene 1): -0.195169
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 8.2
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.6
#> -----------------------------------------------------
#> Iteration 1440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 658.774
#> delta (gene 1): 1.15565
#> phi (cell 1): 0.770987
#> s (cell 1): 0.390202
#> nu (cell 1): 0.421879
#> theta (batch 1): 0.165539
#> betas:    0.6711  -0.1539   0.3204   0.7001   1.4375   0.2675   0.8454   0.2001   0.4825   0.0092  -0.3039   0.2520
#> 
#> sigma: 0.213851
#> lambda (gene 1): 1.1019
#> epsilon (gene 1): 0.399643
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 8.2
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.6
#> -----------------------------------------------------
#> Iteration 1460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 631.907
#> delta (gene 1): 0.590279
#> phi (cell 1): 0.804024
#> s (cell 1): 0.289093
#> nu (cell 1): 0.434444
#> theta (batch 1): 0.300745
#> betas:    1.2324  -0.2214  -0.1098   0.6435   0.3063   0.9103   0.1286   0.5453   0.0134   0.0059  -0.1133   0.1650
#> 
#> sigma: 0.14816
#> lambda (gene 1): 0.500019
#> epsilon (gene 1): -0.423073
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.1
#> LSphi: 8.1
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -1.5
#> -----------------------------------------------------
#> Iteration 1480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 720.814
#> delta (gene 1): 0.811666
#> phi (cell 1): 0.760653
#> s (cell 1): 0.243504
#> nu (cell 1): 0.433343
#> theta (batch 1): 0.589624
#> betas:    1.0537  -0.1774  -0.3520   0.5960   1.2959  -0.0226   0.7714   0.0121   0.7281  -0.8085   0.0119   0.0443
#> 
#> sigma: 0.219353
#> lambda (gene 1): 1.46467
#> epsilon (gene 1): 0.182109
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.1
#> LSphi: 8.1
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -1.5
#> -----------------------------------------------------
#> Iteration 1500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 765.206
#> delta (gene 1): 0.594397
#> phi (cell 1): 0.62509
#> s (cell 1): 0.349992
#> nu (cell 1): 0.431662
#> theta (batch 1): 0.258653
#> betas:    0.7679  -0.1058   0.7015   0.4304   0.8284   0.6290   0.1584   0.9066  -0.6232   0.2317  -0.0399  -0.6330
#> 
#> sigma: 0.139253
#> lambda (gene 1): 0.537802
#> epsilon (gene 1): -0.431221
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 1520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 752.24
#> delta (gene 1): 0.759991
#> phi (cell 1): 0.624734
#> s (cell 1): 0.615933
#> nu (cell 1): 0.429083
#> theta (batch 1): 0.14049
#> betas:    0.5879  -0.0050   1.0159   0.7980   1.0872  -0.0071   0.1610   0.4704   0.0122  -0.6489  -0.1141  -0.4800
#> 
#> sigma: 0.16688
#> lambda (gene 1): 2.42233
#> epsilon (gene 1): -0.109175
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 1540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 672.661
#> delta (gene 1): 0.654832
#> phi (cell 1): 0.673729
#> s (cell 1): 0.739511
#> nu (cell 1): 0.43916
#> theta (batch 1): 0.227748
#> betas:    1.0908  -0.1704   1.2415   0.3378   0.9018   0.4736   0.4675  -0.0291   0.3743  -0.3952   0.0040  -0.0134
#> 
#> sigma: 0.139373
#> lambda (gene 1): 0.267101
#> epsilon (gene 1): -0.2769
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 1560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 840.932
#> delta (gene 1): 0.688165
#> phi (cell 1): 0.733378
#> s (cell 1): 0.551345
#> nu (cell 1): 0.408413
#> theta (batch 1): 0.229235
#> betas:    1.3538  -0.2569   1.0572   0.8099   0.2630   0.6227   0.6961   0.3293  -0.2161   0.4294  -0.0421  -0.4567
#> 
#> sigma: 0.206399
#> lambda (gene 1): 1.84526
#> epsilon (gene 1): -0.277565
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 7.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 1580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 742.626
#> delta (gene 1): 1.10493
#> phi (cell 1): 0.685123
#> s (cell 1): 0.304877
#> nu (cell 1): 0.436138
#> theta (batch 1): 0.245565
#> betas:    0.3773  -0.0052   0.5290   0.7136   0.7198   1.1618  -0.2787   0.7606  -0.2081  -0.4058  -0.0440  -0.5973
#> 
#> sigma: 0.163352
#> lambda (gene 1): 0.234357
#> epsilon (gene 1): 0.327459
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 7.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 1600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 672.79
#> delta (gene 1): 0.619995
#> phi (cell 1): 0.695418
#> s (cell 1): 0.354708
#> nu (cell 1): 0.40913
#> theta (batch 1): 0.139736
#> betas:    0.8592  -0.1693   0.2857   0.9816   0.9831   0.2233   0.7425   0.1908   0.7687  -1.0448   0.5971   0.0311
#> 
#> sigma: 0.212939
#> lambda (gene 1): 1.05188
#> epsilon (gene 1): -0.22639
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 7.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 1620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 629.229
#> delta (gene 1): 0.647719
#> phi (cell 1): 0.674468
#> s (cell 1): 0.379988
#> nu (cell 1): 0.419819
#> theta (batch 1): 0.0995592
#> betas:    0.7974  -0.0414   0.3568   0.9823   0.3124   0.6068   0.6958  -0.3763   0.0218  -0.1422  -0.3798  -0.7201
#> 
#> sigma: 0.156666
#> lambda (gene 1): 1.17981
#> epsilon (gene 1): -0.295767
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 7.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 1640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 740.601
#> delta (gene 1): 0.637581
#> phi (cell 1): 0.657727
#> s (cell 1): 0.260026
#> nu (cell 1): 0.422469
#> theta (batch 1): 0.313503
#> betas:    0.7957  -0.1184   0.3590   0.8952   0.3759   0.9178   0.3235   0.2398   0.2127  -0.0437  -0.2822   0.0688
#> 
#> sigma: 0.189334
#> lambda (gene 1): 0.168359
#> epsilon (gene 1): -0.448651
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 7.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 1660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 704.073
#> delta (gene 1): 0.50707
#> phi (cell 1): 0.622009
#> s (cell 1): 0.813787
#> nu (cell 1): 0.416395
#> theta (batch 1): 0.42891
#> betas:    0.2340   0.0108   0.7606   0.6574   0.3792   1.0641   0.3155   0.6677  -0.1967  -0.3685  -0.3262  -0.0421
#> 
#> sigma: 0.184721
#> lambda (gene 1): 1.4745
#> epsilon (gene 1): -0.387286
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.3
#> LSphi: 7.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 1680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 611.629
#> delta (gene 1): 0.805595
#> phi (cell 1): 0.650249
#> s (cell 1): 0.75667
#> nu (cell 1): 0.420679
#> theta (batch 1): 0.128439
#> betas:   -0.0273   0.0361   0.4475   0.8219   0.7132   0.9829   0.6326   0.1006  -0.0348   0.2230  -0.9783  -0.2336
#> 
#> sigma: 0.165938
#> lambda (gene 1): 3.01984
#> epsilon (gene 1): 0.13644
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.3
#> LSphi: 7.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 1700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 759.8
#> delta (gene 1): 0.622433
#> phi (cell 1): 0.735124
#> s (cell 1): 0.557044
#> nu (cell 1): 0.428804
#> theta (batch 1): 0.328575
#> betas:    0.9756  -0.0430   0.5799   0.3892   0.5948   0.2913   0.8371  -0.2687  -0.3293   0.0183  -0.5006  -0.8147
#> 
#> sigma: 0.213429
#> lambda (gene 1): 0.905107
#> epsilon (gene 1): -0.292827
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7.6
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 1720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 712.724
#> delta (gene 1): 0.576032
#> phi (cell 1): 0.723926
#> s (cell 1): 0.503848
#> nu (cell 1): 0.421915
#> theta (batch 1): 0.152791
#> betas:    0.2009   0.0086   0.9926   1.0558   0.0833   0.9615   0.5475  -0.1144   0.1779  -0.1266  -0.4015  -0.2847
#> 
#> sigma: 0.162599
#> lambda (gene 1): 0.616919
#> epsilon (gene 1): -0.376244
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7.6
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 1740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 802.993
#> delta (gene 1): 0.756028
#> phi (cell 1): 0.737131
#> s (cell 1): 0.78622
#> nu (cell 1): 0.429117
#> theta (batch 1): 0.0907174
#> betas:   -0.0237   0.1066   0.5414   0.3622   1.4482   0.4029   0.4613  -0.0194   0.1167  -0.5428  -0.7543  -0.2973
#> 
#> sigma: 0.156326
#> lambda (gene 1): 0.847854
#> epsilon (gene 1): 0.189931
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7.6
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 1760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 945.675
#> delta (gene 1): 0.584723
#> phi (cell 1): 0.765217
#> s (cell 1): 0.547117
#> nu (cell 1): 0.417447
#> theta (batch 1): 0.228165
#> betas:    1.5956  -0.1996   0.6951   0.6173   0.2074   0.7012   0.0581   0.2676  -0.3303   0.0233  -0.4204  -0.2175
#> 
#> sigma: 0.158333
#> lambda (gene 1): 0.875729
#> epsilon (gene 1): -0.306509
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 7.5
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 1780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 531.061
#> delta (gene 1): 0.727407
#> phi (cell 1): 0.746874
#> s (cell 1): 0.343098
#> nu (cell 1): 0.425018
#> theta (batch 1): 0.264762
#> betas:    0.7825  -0.0122   0.2996   1.0679   0.8094  -0.0157   0.5057  -0.0279   0.1756  -0.1850  -0.6788  -0.6760
#> 
#> sigma: 0.1597
#> lambda (gene 1): 2.41923
#> epsilon (gene 1): -0.227302
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 7.5
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 1800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 647.695
#> delta (gene 1): 0.598974
#> phi (cell 1): 0.759303
#> s (cell 1): 1.44081
#> nu (cell 1): 0.429856
#> theta (batch 1): 0.381063
#> betas:    0.5552  -0.0681   0.3037   0.6349   1.2781   0.3862   0.2674   0.2231   0.0792  -0.1725  -0.4439   0.0128
#> 
#> sigma: 0.173945
#> lambda (gene 1): 0.549643
#> epsilon (gene 1): -0.281575
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7.4
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 1820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 649.778
#> delta (gene 1): 0.59596
#> phi (cell 1): 0.744051
#> s (cell 1): 0.427317
#> nu (cell 1): 0.422018
#> theta (batch 1): 0.147812
#> betas:    0.9404  -0.0935   0.1344   0.7990   0.7676   0.4819   0.5563   0.0237   0.0490  -0.0071  -0.6404  -0.1099
#> 
#> sigma: 0.130244
#> lambda (gene 1): 1.30198
#> epsilon (gene 1): -0.448231
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7.4
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 1840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 829.756
#> delta (gene 1): 0.602195
#> phi (cell 1): 0.746001
#> s (cell 1): 0.642623
#> nu (cell 1): 0.427305
#> theta (batch 1): 0.388523
#> betas:    0.4159  -0.0454   0.2717   1.0212   1.1394   0.6989   0.4308  -0.0916   0.4598  -0.1255  -0.3264  -0.3751
#> 
#> sigma: 0.141421
#> lambda (gene 1): 0.31791
#> epsilon (gene 1): -0.521516
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7.4
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 1860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 750.964
#> delta (gene 1): 0.935934
#> phi (cell 1): 0.673258
#> s (cell 1): 0.515231
#> nu (cell 1): 0.422086
#> theta (batch 1): 0.265189
#> betas:    0.5980  -0.0353   0.4138   0.9473   1.0796  -0.0349   0.2669   0.3153   0.4180  -0.6805  -0.0424  -0.4786
#> 
#> sigma: 0.178286
#> lambda (gene 1): 0.260751
#> epsilon (gene 1): 0.0261946
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 7.3
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 1880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 671.792
#> delta (gene 1): 0.666414
#> phi (cell 1): 0.712222
#> s (cell 1): 0.380761
#> nu (cell 1): 0.420967
#> theta (batch 1): 0.254424
#> betas:    0.9223  -0.0644   0.3828   1.0243   0.6599   0.6597  -0.0504  -0.1053   0.5186  -0.7820   0.1834  -0.7073
#> 
#> sigma: 0.12415
#> lambda (gene 1): 1.46701
#> epsilon (gene 1): -0.41261
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 7.3
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 1900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 722.183
#> delta (gene 1): 0.684336
#> phi (cell 1): 0.763982
#> s (cell 1): 0.51989
#> nu (cell 1): 0.43608
#> theta (batch 1): 0.0778678
#> betas:    1.3397  -0.1619   0.2268   0.8363   0.1340   0.7075   0.1899   0.5189  -0.2714  -0.1667  -0.3234   0.0423
#> 
#> sigma: 0.165004
#> lambda (gene 1): 0.713746
#> epsilon (gene 1): -0.173625
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7.2
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 1920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 604.985
#> delta (gene 1): 0.454133
#> phi (cell 1): 0.70403
#> s (cell 1): 0.446702
#> nu (cell 1): 0.42211
#> theta (batch 1): 0.129673
#> betas:    0.2626  -0.1570   1.1581  -0.1901   0.6994   1.5453   0.4590   0.3270   0.7451   0.2992  -0.6396   0.6724
#> 
#> sigma: 0.165253
#> lambda (gene 1): 1.14547
#> epsilon (gene 1): -0.506586
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7.2
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 1940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 688.095
#> delta (gene 1): 0.684054
#> phi (cell 1): 0.793785
#> s (cell 1): 0.537637
#> nu (cell 1): 0.42105
#> theta (batch 1): 0.185443
#> betas:    1.2108  -0.1313   0.4768   0.0970   1.2205  -0.2065   0.7647   0.3496  -0.3888   0.1166  -0.7426   0.0746
#> 
#> sigma: 0.171825
#> lambda (gene 1): 1.38433
#> epsilon (gene 1): -0.102125
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7.2
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 1960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 752.675
#> delta (gene 1): 0.645534
#> phi (cell 1): 0.738499
#> s (cell 1): 0.60903
#> nu (cell 1): 0.409752
#> theta (batch 1): 0.0522116
#> betas:    0.3372  -0.0622  -0.4475   0.5279   0.7525   1.5377   0.4589   0.4728  -0.0265   0.1815  -0.2933  -0.2630
#> 
#> sigma: 0.127475
#> lambda (gene 1): 1.46241
#> epsilon (gene 1): -0.264724
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 7.1
#> LSnu (cell 1): -6.41822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 1980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 757.127
#> delta (gene 1): 0.74548
#> phi (cell 1): 0.648838
#> s (cell 1): 0.392754
#> nu (cell 1): 0.421109
#> theta (batch 1): 0.0618186
#> betas:    0.6157  -0.0741   0.1605   0.6585   0.6665   1.2131   0.4361   0.4193  -0.1309   0.1222  -0.2803  -0.3689
#> 
#> sigma: 0.167583
#> lambda (gene 1): 2.8698
#> epsilon (gene 1): -0.171628
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 7.1
#> LSnu (cell 1): -6.41822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 2000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 701.553
#> delta (gene 1): 0.998727
#> phi (cell 1): 0.632522
#> s (cell 1): 0.605726
#> nu (cell 1): 0.420683
#> theta (batch 1): 0.0987985
#> betas:    0.3646  -0.1481   0.4000   0.2639   1.0163   1.2044   0.3173   0.8409   0.2272  -0.1073   0.1958   0.2589
#> 
#> sigma: 0.168146
#> lambda (gene 1): 0.631481
#> epsilon (gene 1): 0.180488
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -6.51822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 831.276
#> delta (gene 1): 0.567679
#> phi (cell 1): 0.818522
#> s (cell 1): 0.613376
#> nu (cell 1): 0.415182
#> theta (batch 1): 0.267135
#> betas:    1.1860  -0.2065   0.2394   0.2392   0.8066   0.5532   0.6637   0.5518  -0.1651   0.0403  -0.3977   0.3544
#> 
#> sigma: 0.149999
#> lambda (gene 1): 1.25243
#> epsilon (gene 1): -0.243586
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -6.51822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 960.311
#> delta (gene 1): 0.504644
#> phi (cell 1): 0.973931
#> s (cell 1): 0.566536
#> nu (cell 1): 0.413067
#> theta (batch 1): 0.433263
#> betas:    0.6404   0.0308   0.9728   0.7012   0.8368   0.3603   0.8569  -0.0521  -0.1537  -0.2352  -0.3647  -0.8071
#> 
#> sigma: 0.178899
#> lambda (gene 1): 1.15611
#> epsilon (gene 1): -0.617089
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -6.51822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 662.177
#> delta (gene 1): 0.683584
#> phi (cell 1): 0.876432
#> s (cell 1): 0.357666
#> nu (cell 1): 0.419623
#> theta (batch 1): 0.378883
#> betas:    1.2893  -0.1296   0.7714  -0.0943   0.4735   0.8163   0.3264  -0.0389   0.1956  -0.0356  -0.5422  -0.3377
#> 
#> sigma: 0.141836
#> lambda (gene 1): 0.996907
#> epsilon (gene 1): -0.171027
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.41822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 2080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 621.897
#> delta (gene 1): 0.956426
#> phi (cell 1): 0.813453
#> s (cell 1): 0.693033
#> nu (cell 1): 0.430374
#> theta (batch 1): 0.332888
#> betas:    0.9696  -0.1848   0.8188   0.1382   0.9513   0.8114   0.8461   0.3966   0.2185  -0.3357   0.0633   0.0334
#> 
#> sigma: 0.138902
#> lambda (gene 1): 0.459882
#> epsilon (gene 1): 0.307814
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.41822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 2100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 632.419
#> delta (gene 1): 0.948901
#> phi (cell 1): 0.8108
#> s (cell 1): 0.5299
#> nu (cell 1): 0.414029
#> theta (batch 1): 0.239393
#> betas:    1.5223  -0.2277  -0.1849   0.1564   0.6982   0.3903   0.9104  -0.1692   0.2179  -0.0182  -0.2135  -0.1378
#> 
#> sigma: 0.239619
#> lambda (gene 1): 5.56161
#> epsilon (gene 1): 0.0834028
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.6
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 547.544
#> delta (gene 1): 0.540259
#> phi (cell 1): 0.72956
#> s (cell 1): 1.05165
#> nu (cell 1): 0.417442
#> theta (batch 1): 0.38631
#> betas:    0.7349  -0.0450   0.2812   0.0176   1.0542   0.4587   0.5680  -0.0319   0.3471  -0.5735  -0.4711  -0.1399
#> 
#> sigma: 0.131714
#> lambda (gene 1): 1.81513
#> epsilon (gene 1): -0.235348
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.6
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 772.771
#> delta (gene 1): 0.586361
#> phi (cell 1): 0.781454
#> s (cell 1): 0.369868
#> nu (cell 1): 0.428318
#> theta (batch 1): 0.500464
#> betas:    1.0981  -0.1275   0.3128   1.0098   0.8325   0.9630  -0.0690   0.4849  -0.2184  -0.1760   0.1936  -0.8092
#> 
#> sigma: 0.158197
#> lambda (gene 1): 1.1172
#> epsilon (gene 1): -0.425105
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.6
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -6.31822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 646.012
#> delta (gene 1): 0.540029
#> phi (cell 1): 0.796251
#> s (cell 1): 0.446617
#> nu (cell 1): 0.414636
#> theta (batch 1): 0.347264
#> betas:    0.9186  -0.1833   0.9416   0.5725   0.3462   0.9210   0.3133   0.7897   0.2895  -0.3208   0.0785  -0.0375
#> 
#> sigma: 0.144119
#> lambda (gene 1): 0.767074
#> epsilon (gene 1): -0.306266
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.5
#> -----------------------------------------------------
#> Iteration 2180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 838.673
#> delta (gene 1): 0.633181
#> phi (cell 1): 0.8996
#> s (cell 1): 0.44677
#> nu (cell 1): 0.417625
#> theta (batch 1): 0.253241
#> betas:    0.9102  -0.1275   0.0149   0.8534   0.6893   0.8204   0.3594   0.0593   0.3585  -0.5646   0.1323  -0.4134
#> 
#> sigma: 0.195934
#> lambda (gene 1): 0.283481
#> epsilon (gene 1): -0.114564
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.5
#> -----------------------------------------------------
#> Iteration 2200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 687.859
#> delta (gene 1): 0.916601
#> phi (cell 1): 0.842146
#> s (cell 1): 0.497657
#> nu (cell 1): 0.417586
#> theta (batch 1): 0.289632
#> betas:    0.4194  -0.0513   0.4112   0.8221   1.0876   0.4586   0.8106   0.2672   0.1521  -0.2230  -0.1244  -0.3115
#> 
#> sigma: 0.139818
#> lambda (gene 1): 0.891362
#> epsilon (gene 1): 0.191285
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 508.591
#> delta (gene 1): 0.85135
#> phi (cell 1): 0.787188
#> s (cell 1): 1.48081
#> nu (cell 1): 0.424092
#> theta (batch 1): 0.273143
#> betas:    0.7720  -0.0701   0.7044   1.6100   0.1583   0.4916   0.9153  -0.0992   0.4739  -0.7417   0.1609  -0.3759
#> 
#> sigma: 0.167756
#> lambda (gene 1): 0.536396
#> epsilon (gene 1): 0.0274062
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 811.573
#> delta (gene 1): 0.71085
#> phi (cell 1): 0.715645
#> s (cell 1): 0.412681
#> nu (cell 1): 0.438982
#> theta (batch 1): 0.216462
#> betas:    1.1863  -0.0316   0.2357   0.9190   0.5058   0.5508   0.0117   0.1931  -0.1631   0.0032  -0.9674  -0.6670
#> 
#> sigma: 0.183431
#> lambda (gene 1): 0.621767
#> epsilon (gene 1): -0.067719
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 738.019
#> delta (gene 1): 0.609241
#> phi (cell 1): 0.640411
#> s (cell 1): 0.770243
#> nu (cell 1): 0.418929
#> theta (batch 1): 0.279361
#> betas:    1.7000  -0.2336   0.1250   0.4214   0.3784   0.7634   0.0730   0.2022   0.1122  -0.1081  -0.3750   0.1719
#> 
#> sigma: 0.168704
#> lambda (gene 1): 0.873325
#> epsilon (gene 1): -0.363461
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 2280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 660.086
#> delta (gene 1): 0.666398
#> phi (cell 1): 0.719506
#> s (cell 1): 0.41701
#> nu (cell 1): 0.447328
#> theta (batch 1): 0.301985
#> betas:    0.7246  -0.1249   0.7557   0.2106   0.6828   0.9456   0.5702   0.0076   0.3071   0.0649  -0.2856   0.0373
#> 
#> sigma: 0.123727
#> lambda (gene 1): 1.08661
#> epsilon (gene 1): -0.291981
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 2300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 645.378
#> delta (gene 1): 0.849889
#> phi (cell 1): 0.759633
#> s (cell 1): 0.665975
#> nu (cell 1): 0.44947
#> theta (batch 1): 0.264676
#> betas:    1.0171  -0.0952   0.4753   0.5363   0.7097   0.7025   0.5191   0.2479  -0.4120   0.1396   0.0005  -0.7302
#> 
#> sigma: 0.152932
#> lambda (gene 1): 0.40235
#> epsilon (gene 1): -0.256104
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 570.939
#> delta (gene 1): 0.969342
#> phi (cell 1): 0.795594
#> s (cell 1): 0.304032
#> nu (cell 1): 0.426632
#> theta (batch 1): 0.181702
#> betas:    1.3266  -0.0677   0.3404  -0.0836   0.4497   0.6339   0.0635   0.1366   0.1980  -0.4793  -0.5711  -0.6132
#> 
#> sigma: 0.213126
#> lambda (gene 1): 0.948957
#> epsilon (gene 1): 0.115696
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 727.651
#> delta (gene 1): 0.636333
#> phi (cell 1): 0.9567
#> s (cell 1): 0.44052
#> nu (cell 1): 0.422599
#> theta (batch 1): 0.212931
#> betas:    0.8117  -0.2128   0.4243   0.7391   1.1595   0.2654   0.9370   0.4309  -0.0714   0.4887  -0.1927   0.0626
#> 
#> sigma: 0.165699
#> lambda (gene 1): 0.703813
#> epsilon (gene 1): -0.272405
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 2360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 753.739
#> delta (gene 1): 0.720275
#> phi (cell 1): 0.965182
#> s (cell 1): 0.582766
#> nu (cell 1): 0.420267
#> theta (batch 1): 0.147069
#> betas:    0.6706  -0.0709   0.3027   0.6444   0.6748   1.0535   0.3318   0.1643   0.1338  -0.0894  -0.1595  -0.4621
#> 
#> sigma: 0.176755
#> lambda (gene 1): 0.740472
#> epsilon (gene 1): -0.249117
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 2380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 831.454
#> delta (gene 1): 0.884712
#> phi (cell 1): 0.967225
#> s (cell 1): 0.518791
#> nu (cell 1): 0.413712
#> theta (batch 1): 0.061764
#> betas:    1.1944  -0.1393   0.2264   0.3257   1.1162   0.5958   0.0524   0.2719   0.0803  -0.3920  -0.2225   0.0709
#> 
#> sigma: 0.143545
#> lambda (gene 1): 1.19024
#> epsilon (gene 1): 0.0540676
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 2400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 621.868
#> delta (gene 1): 0.793508
#> phi (cell 1): 0.998704
#> s (cell 1): 0.430838
#> nu (cell 1): 0.423427
#> theta (batch 1): 0.0580407
#> betas:    0.6824   0.0307  -0.2171   0.8856   0.5567   0.8593   0.4634  -0.3759   0.4035  -0.6645  -0.5330  -0.4193
#> 
#> sigma: 0.200031
#> lambda (gene 1): 1.50763
#> epsilon (gene 1): -0.305559
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 2420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 779.447
#> delta (gene 1): 0.865957
#> phi (cell 1): 0.886216
#> s (cell 1): 0.538544
#> nu (cell 1): 0.418943
#> theta (batch 1): 0.0333373
#> betas:    1.0407  -0.1075   0.9555   0.8620  -0.2545   1.0806   0.4779  -0.2677   0.3806  -0.0253  -0.7823   0.1277
#> 
#> sigma: 0.186
#> lambda (gene 1): 2.99331
#> epsilon (gene 1): -0.323891
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 2440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1014.92
#> delta (gene 1): 1.24768
#> phi (cell 1): 0.887238
#> s (cell 1): 0.653991
#> nu (cell 1): 0.434952
#> theta (batch 1): 0.11486
#> betas:    0.7988   0.0644   0.4856   0.6081   0.3014   0.4519   0.7184  -0.3011  -0.0128  -0.8007  -0.5304  -0.5476
#> 
#> sigma: 0.197667
#> lambda (gene 1): 1.33612
#> epsilon (gene 1): 0.277915
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 2460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 719.996
#> delta (gene 1): 0.670511
#> phi (cell 1): 0.814758
#> s (cell 1): 0.526673
#> nu (cell 1): 0.412102
#> theta (batch 1): 0.215348
#> betas:    1.5823  -0.2682   0.2693  -0.3151   0.6031   0.7406  -0.0095   0.2436   0.3334  -0.2074   0.2666  -0.1754
#> 
#> sigma: 0.162731
#> lambda (gene 1): 0.30142
#> epsilon (gene 1): -0.371116
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 2480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 621.242
#> delta (gene 1): 0.786384
#> phi (cell 1): 0.892032
#> s (cell 1): 0.862842
#> nu (cell 1): 0.419769
#> theta (batch 1): 0.258001
#> betas:    0.5358  -0.0898   0.7781   0.9007   0.7131   0.7846   0.4263   0.2539   0.1312   0.0670  -0.3583  -0.0884
#> 
#> sigma: 0.146726
#> lambda (gene 1): 2.19899
#> epsilon (gene 1): 0.00173689
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 2500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 782.813
#> delta (gene 1): 0.718578
#> phi (cell 1): 0.800714
#> s (cell 1): 0.828246
#> nu (cell 1): 0.416869
#> theta (batch 1): 0.209973
#> betas:    1.2054   0.0233  -0.3409   0.7312   0.1641   1.2744  -0.5003   0.4659  -0.7157   0.2730  -1.3398  -0.7768
#> 
#> sigma: 0.196758
#> lambda (gene 1): 0.591951
#> epsilon (gene 1): -0.00975186
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 2520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 763.881
#> delta (gene 1): 0.620685
#> phi (cell 1): 0.746921
#> s (cell 1): 0.253205
#> nu (cell 1): 0.424089
#> theta (batch 1): 0.421804
#> betas:    0.8585  -0.0243   1.4785   0.5556   0.5468  -0.0089   1.1927  -0.5458   0.2273  -0.1591  -0.7216  -0.6527
#> 
#> sigma: 0.160671
#> lambda (gene 1): 1.23081
#> epsilon (gene 1): -0.0642045
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 2540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 715.119
#> delta (gene 1): 1.00564
#> phi (cell 1): 0.795788
#> s (cell 1): 0.343154
#> nu (cell 1): 0.417369
#> theta (batch 1): 0.0887734
#> betas:    1.4106  -0.1186   0.5676   0.2108   1.1535   0.3364  -0.1024   0.2709  -0.0622   0.3774  -1.0936  -0.2434
#> 
#> sigma: 0.20527
#> lambda (gene 1): 0.874126
#> epsilon (gene 1): 0.230216
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 2560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 662.023
#> delta (gene 1): 0.525451
#> phi (cell 1): 0.832816
#> s (cell 1): 0.545927
#> nu (cell 1): 0.411685
#> theta (batch 1): 0.104993
#> betas:    0.3414   0.0049   1.0754   0.6597   1.0066   0.4221   0.6651  -0.0040   0.1473  -0.0087  -0.2448  -0.7500
#> 
#> sigma: 0.151063
#> lambda (gene 1): 0.721423
#> epsilon (gene 1): -0.649169
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 2580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.167
#> delta (gene 1): 0.834415
#> phi (cell 1): 0.792398
#> s (cell 1): 0.371772
#> nu (cell 1): 0.429632
#> theta (batch 1): 0.179065
#> betas:    1.1167  -0.1333   0.2884   0.4898   0.9793   0.5575   0.4091  -0.1523   0.5451  -0.3057  -0.1119  -0.2981
#> 
#> sigma: 0.174905
#> lambda (gene 1): 1.06597
#> epsilon (gene 1): -0.233679
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 2600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 624.049
#> delta (gene 1): 0.510827
#> phi (cell 1): 0.837935
#> s (cell 1): 0.472402
#> nu (cell 1): 0.422835
#> theta (batch 1): 0.0457781
#> betas:    0.5407  -0.0278   0.0823   1.0188   0.2013   0.9332   0.8505  -0.2020   0.2665  -0.4219  -0.3170  -0.1374
#> 
#> sigma: 0.174328
#> lambda (gene 1): 1.07699
#> epsilon (gene 1): -0.522237
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 2620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 624.106
#> delta (gene 1): 0.696463
#> phi (cell 1): 0.795403
#> s (cell 1): 0.602482
#> nu (cell 1): 0.424287
#> theta (batch 1): 0.0309949
#> betas:    0.9635  -0.0560   0.6420   0.1389   0.8227   0.4757   0.3353   0.0971   0.0492  -0.3355  -0.4898  -0.4202
#> 
#> sigma: 0.184811
#> lambda (gene 1): 1.24316
#> epsilon (gene 1): -0.240861
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 2640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 846.581
#> delta (gene 1): 0.89255
#> phi (cell 1): 0.915716
#> s (cell 1): 0.445626
#> nu (cell 1): 0.419536
#> theta (batch 1): 0.070808
#> betas:    0.7993  -0.0786   0.3550   1.0367   0.7669   0.1637   1.0823   0.0268  -0.3754   0.2453  -0.7138  -0.1535
#> 
#> sigma: 0.159919
#> lambda (gene 1): 1.51854
#> epsilon (gene 1): 0.109716
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 2660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 465.617
#> delta (gene 1): 0.943887
#> phi (cell 1): 0.90021
#> s (cell 1): 0.445009
#> nu (cell 1): 0.418994
#> theta (batch 1): 0.0151478
#> betas:   -0.2065   0.2048   1.2167   1.0681   0.0541   0.5910   0.2117  -0.1486  -0.4900  -0.5284  -0.4595  -0.9041
#> 
#> sigma: 0.176947
#> lambda (gene 1): 0.476217
#> epsilon (gene 1): 0.131002
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.5
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.7
#> -----------------------------------------------------
#> Iteration 2680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 691.373
#> delta (gene 1): 0.700782
#> phi (cell 1): 0.761593
#> s (cell 1): 0.607472
#> nu (cell 1): 0.423674
#> theta (batch 1): 0.121265
#> betas:    0.7905  -0.0516   0.6653   1.0422   0.4291  -0.1521   0.5740   0.0274  -0.1390  -0.0583  -0.9448   0.3644
#> 
#> sigma: 0.156057
#> lambda (gene 1): 1.29918
#> epsilon (gene 1): -0.298907
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.5
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.7
#> -----------------------------------------------------
#> Iteration 2700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 711.72
#> delta (gene 1): 0.934163
#> phi (cell 1): 0.635977
#> s (cell 1): 0.463621
#> nu (cell 1): 0.413342
#> theta (batch 1): 0.320513
#> betas:    0.5174   0.0102   0.8245   0.4903   0.9085   0.3749   0.4797   0.1526   0.0864  -0.2923  -0.4969  -0.7449
#> 
#> sigma: 0.210703
#> lambda (gene 1): 1.27403
#> epsilon (gene 1): 0.04032
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 2720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 797.406
#> delta (gene 1): 0.771811
#> phi (cell 1): 0.665287
#> s (cell 1): 0.60377
#> nu (cell 1): 0.413394
#> theta (batch 1): 0.296801
#> betas:    0.2573  -0.0187   0.6559   0.2080   0.7916   1.3471   0.3594  -0.0448   0.3714  -0.2866  -0.0823  -0.6865
#> 
#> sigma: 0.257276
#> lambda (gene 1): 1.60736
#> epsilon (gene 1): -0.167084
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 2740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 646.724
#> delta (gene 1): 0.571596
#> phi (cell 1): 0.732891
#> s (cell 1): 0.850768
#> nu (cell 1): 0.425391
#> theta (batch 1): 0.109293
#> betas:    0.9204  -0.1397   0.1895   0.5817   0.3887   0.6406   0.5972   0.2853   0.1842  -0.1151  -0.3017  -0.1312
#> 
#> sigma: 0.15512
#> lambda (gene 1): 1.80439
#> epsilon (gene 1): -0.474918
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 2760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 705.134
#> delta (gene 1): 1.14138
#> phi (cell 1): 0.787281
#> s (cell 1): 0.484111
#> nu (cell 1): 0.414741
#> theta (batch 1): 0.154984
#> betas:    0.1140   0.1782   1.1957   0.5115   1.0453   0.4066   0.1946   0.4853  -0.3793  -0.3033  -0.8122  -0.9043
#> 
#> sigma: 0.152966
#> lambda (gene 1): 1.1501
#> epsilon (gene 1): 0.293946
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 2780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 592.25
#> delta (gene 1): 0.578761
#> phi (cell 1): 0.667012
#> s (cell 1): 0.792359
#> nu (cell 1): 0.429083
#> theta (batch 1): 0.30791
#> betas:    0.7261   0.0164   0.4017   0.4604   0.5947   0.6842   0.3777   0.4874  -0.8978   0.2881  -0.8865  -0.4404
#> 
#> sigma: 0.186845
#> lambda (gene 1): 0.828371
#> epsilon (gene 1): -0.329311
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.1
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 2800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 702.7
#> delta (gene 1): 0.51441
#> phi (cell 1): 0.781232
#> s (cell 1): 0.633613
#> nu (cell 1): 0.420273
#> theta (batch 1): 0.301488
#> betas:    0.3346   0.1374   0.4348   0.0962   0.7079   0.6985   0.6280  -0.0929   0.1920  -0.5714  -0.7378  -0.9019
#> 
#> sigma: 0.148647
#> lambda (gene 1): 1.08695
#> epsilon (gene 1): -0.358138
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.2
#> LSphi: 6.6
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 2820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 626.426
#> delta (gene 1): 1.00164
#> phi (cell 1): 0.81343
#> s (cell 1): 0.42759
#> nu (cell 1): 0.42169
#> theta (batch 1): 0.258283
#> betas:    1.5592  -0.1561   0.1553   1.1051   0.6221   0.4130   0.2068   0.5669  -0.4742   0.2411  -0.6415  -0.3701
#> 
#> sigma: 0.137314
#> lambda (gene 1): 1.02655
#> epsilon (gene 1): 0.177754
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.2
#> LSphi: 6.6
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 2840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.793
#> delta (gene 1): 0.545878
#> phi (cell 1): 0.703321
#> s (cell 1): 0.526354
#> nu (cell 1): 0.430068
#> theta (batch 1): 0.113894
#> betas:    0.8161  -0.0405   0.4257   1.3121   0.1139   0.9007   0.4964  -0.2216   0.3864   0.3187  -1.0794  -0.1205
#> 
#> sigma: 0.19699
#> lambda (gene 1): 0.777038
#> epsilon (gene 1): -0.531465
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.2
#> LSphi: 6.6
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 2860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 848.251
#> delta (gene 1): 0.797412
#> phi (cell 1): 0.682874
#> s (cell 1): 0.259012
#> nu (cell 1): 0.406491
#> theta (batch 1): 0.0474569
#> betas:    0.7621  -0.0191   0.5741  -0.1600   1.1386   0.6249   0.1052   0.3878  -0.0970  -0.3914  -0.4368  -0.5151
#> 
#> sigma: 0.151689
#> lambda (gene 1): 1.45087
#> epsilon (gene 1): -0.091433
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 2880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 759.479
#> delta (gene 1): 0.829324
#> phi (cell 1): 0.579978
#> s (cell 1): 0.480027
#> nu (cell 1): 0.424579
#> theta (batch 1): 0.090351
#> betas:    0.6288  -0.0165   0.5311   1.3185   0.4421   0.8238   0.3437   0.0206   0.0048  -0.4302  -0.1117  -0.6236
#> 
#> sigma: 0.154051
#> lambda (gene 1): 1.63976
#> epsilon (gene 1): -0.0844617
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 2900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 613.054
#> delta (gene 1): 0.700363
#> phi (cell 1): 0.546248
#> s (cell 1): 0.67792
#> nu (cell 1): 0.419824
#> theta (batch 1): 0.248861
#> betas:    1.4626  -0.1415   0.4896  -0.6622   0.8470   0.4746   0.2474   0.0437   0.0668   0.0018  -0.8281  -0.2569
#> 
#> sigma: 0.190403
#> lambda (gene 1): 0.583495
#> epsilon (gene 1): -0.493944
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 2920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 553.175
#> delta (gene 1): 1.36972
#> phi (cell 1): 0.625031
#> s (cell 1): 0.854132
#> nu (cell 1): 0.412002
#> theta (batch 1): 0.206436
#> betas:    0.7146  -0.0750   0.2025   1.2757   0.4442   0.4986   0.6747   0.5008  -0.2815   0.2557  -0.6081   0.0151
#> 
#> sigma: 0.176116
#> lambda (gene 1): 0.875574
#> epsilon (gene 1): 0.186372
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 2940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 578.087
#> delta (gene 1): 0.661573
#> phi (cell 1): 0.64726
#> s (cell 1): 0.337654
#> nu (cell 1): 0.410454
#> theta (batch 1): 0.0510555
#> betas:    1.1494  -0.0749  -0.0666   0.8724   0.3802   0.7367   0.1729   0.0929   0.1209  -0.2287  -0.9217   0.0183
#> 
#> sigma: 0.117866
#> lambda (gene 1): 0.693551
#> epsilon (gene 1): -0.468629
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 2960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 759.522
#> delta (gene 1): 0.600509
#> phi (cell 1): 0.723089
#> s (cell 1): 0.308195
#> nu (cell 1): 0.421273
#> theta (batch 1): 0.107661
#> betas:    0.6934  -0.0004   0.6936   0.7502   0.7140   0.2103   0.7454  -0.0637  -0.0994  -0.1955  -0.6599  -0.2262
#> 
#> sigma: 0.208018
#> lambda (gene 1): 0.157551
#> epsilon (gene 1): -0.336782
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 2980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 810.665
#> delta (gene 1): 0.740428
#> phi (cell 1): 0.764147
#> s (cell 1): 0.394377
#> nu (cell 1): 0.429121
#> theta (batch 1): 0.0862311
#> betas:    0.4355   0.0418   0.9629   0.5398   1.0518   0.1628   1.0051  -0.1569  -0.0502   0.0026  -0.6813  -0.5644
#> 
#> sigma: 0.177737
#> lambda (gene 1): 0.17829
#> epsilon (gene 1): -0.222653
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 3000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 739.863
#> delta (gene 1): 0.531501
#> phi (cell 1): 0.789426
#> s (cell 1): 0.376048
#> nu (cell 1): 0.422836
#> theta (batch 1): 0.0884607
#> betas:    1.2497  -0.2231   0.5726   0.4746   0.5432   0.2877   0.7541   0.8676  -0.2404  -0.1442  -0.0512   0.2200
#> 
#> sigma: 0.23767
#> lambda (gene 1): 1.7302
#> epsilon (gene 1): -0.336411
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 736.373
#> delta (gene 1): 0.571528
#> phi (cell 1): 0.740767
#> s (cell 1): 0.396314
#> nu (cell 1): 0.418139
#> theta (batch 1): 0.0937603
#> betas:    0.9709  -0.0566   1.0777  -0.0013   0.5696   0.5187   0.6593   0.4515  -0.3365  -0.0817  -0.4324  -0.5638
#> 
#> sigma: 0.167389
#> lambda (gene 1): 0.832275
#> epsilon (gene 1): -0.335997
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 836.08
#> delta (gene 1): 1.00744
#> phi (cell 1): 0.784914
#> s (cell 1): 0.298406
#> nu (cell 1): 0.421794
#> theta (batch 1): 0.0796855
#> betas:    1.2980  -0.2125  -0.4238   0.6997   0.7264   0.7707   0.5240   0.2721  -0.0238   0.0890   0.0124  -0.1400
#> 
#> sigma: 0.165364
#> lambda (gene 1): 1.30886
#> epsilon (gene 1): 0.0495466
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 618.7
#> delta (gene 1): 0.958093
#> phi (cell 1): 0.758078
#> s (cell 1): 1.22857
#> nu (cell 1): 0.423933
#> theta (batch 1): 0.383076
#> betas:    1.1201  -0.1115   0.2860   1.5915   0.0360   0.5529   0.5595   0.0185  -0.0719  -0.1553  -0.3772  -0.2132
#> 
#> sigma: 0.221677
#> lambda (gene 1): 2.60774
#> epsilon (gene 1): 0.148217
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 7.1
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 3080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 664.645
#> delta (gene 1): 0.641783
#> phi (cell 1): 0.868637
#> s (cell 1): 0.23232
#> nu (cell 1): 0.419697
#> theta (batch 1): 0.0987349
#> betas:    1.0384  -0.1734   0.8624   0.6684   0.7319   0.8971   0.1734   0.4243   0.0240   0.3015  -0.6081   0.1877
#> 
#> sigma: 0.148804
#> lambda (gene 1): 1.1906
#> epsilon (gene 1): -0.26373
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 7.1
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 3100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 642.678
#> delta (gene 1): 0.526146
#> phi (cell 1): 0.708982
#> s (cell 1): 0.507578
#> nu (cell 1): 0.410711
#> theta (batch 1): 0.206057
#> betas:    1.0414  -0.0681   0.8005   0.2446   0.7714   0.6073   0.5971  -0.2807   0.0517  -0.4095  -0.0401  -0.9347
#> 
#> sigma: 0.236526
#> lambda (gene 1): 1.76396
#> epsilon (gene 1): -0.450477
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 714.42
#> delta (gene 1): 0.48729
#> phi (cell 1): 0.700677
#> s (cell 1): 0.6359
#> nu (cell 1): 0.432089
#> theta (batch 1): 0.282072
#> betas:    0.4310  -0.0533   0.8342   0.4442   0.8817   0.8711   0.3804   0.2714  -0.0455   0.5718  -0.7621  -0.2560
#> 
#> sigma: 0.171569
#> lambda (gene 1): 0.744343
#> epsilon (gene 1): -0.616392
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1021.86
#> delta (gene 1): 0.932035
#> phi (cell 1): 0.796441
#> s (cell 1): 0.353919
#> nu (cell 1): 0.429696
#> theta (batch 1): 0.219109
#> betas:    0.5679   0.1105   0.3565   0.7303   0.4884   0.5398   0.6418  -0.5164   0.1728  -0.4301  -0.7486  -1.0400
#> 
#> sigma: 0.209398
#> lambda (gene 1): 0.362872
#> epsilon (gene 1): 0.0845703
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 995.558
#> delta (gene 1): 0.845585
#> phi (cell 1): 0.871706
#> s (cell 1): 0.274637
#> nu (cell 1): 0.411408
#> theta (batch 1): 0.114028
#> betas:    1.1952  -0.2052  -0.2109   0.5061   0.3996   0.5835   0.9097   0.4559  -0.3435   0.5067  -0.4192   0.2057
#> 
#> sigma: 0.114002
#> lambda (gene 1): 1.56023
#> epsilon (gene 1): -0.097466
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.9
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 3180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 800.105
#> delta (gene 1): 0.764793
#> phi (cell 1): 0.841873
#> s (cell 1): 0.405384
#> nu (cell 1): 0.419646
#> theta (batch 1): 0.0446547
#> betas:    1.4728  -0.1832   0.2509   0.3990   0.3745   0.6374   0.0733   0.6071  -0.2162  -0.0523  -0.4800   0.1216
#> 
#> sigma: 0.176479
#> lambda (gene 1): 0.26047
#> epsilon (gene 1): -0.0814497
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.9
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 3200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 785.561
#> delta (gene 1): 0.750465
#> phi (cell 1): 0.864773
#> s (cell 1): 0.364804
#> nu (cell 1): 0.408969
#> theta (batch 1): 0.280465
#> betas:    0.6967  -0.0659   0.3306   0.7571   0.4762   0.9434   0.2206   0.2112  -0.1841   0.3090  -0.8001  -0.2914
#> 
#> sigma: 0.162268
#> lambda (gene 1): 0.703852
#> epsilon (gene 1): -0.0329256
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 699.588
#> delta (gene 1): 0.734791
#> phi (cell 1): 0.799853
#> s (cell 1): 0.282336
#> nu (cell 1): 0.424919
#> theta (batch 1): 0.236861
#> betas:    0.8338  -0.0856   0.3966   0.8727   0.9646   0.1834   0.4125   0.7183  -0.6582  -0.1887  -0.1614  -0.2035
#> 
#> sigma: 0.141475
#> lambda (gene 1): 2.31858
#> epsilon (gene 1): 0.0639715
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 725.757
#> delta (gene 1): 0.576161
#> phi (cell 1): 0.750166
#> s (cell 1): 0.383596
#> nu (cell 1): 0.40262
#> theta (batch 1): 0.0577037
#> betas:    0.8849   0.0317  -0.4821   1.1088   0.6121   0.4593   0.1968   0.3706  -0.3911  -0.2377  -0.9877  -0.2276
#> 
#> sigma: 0.152591
#> lambda (gene 1): 0.891836
#> epsilon (gene 1): -0.373311
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 726.619
#> delta (gene 1): 0.459613
#> phi (cell 1): 0.731598
#> s (cell 1): 1.6334
#> nu (cell 1): 0.432453
#> theta (batch 1): 0.255885
#> betas:    0.9734  -0.1918   0.4684   0.3957   1.0279   0.9643   0.4244   0.3873  -0.2202   0.3562  -0.2500   0.0656
#> 
#> sigma: 0.167524
#> lambda (gene 1): 0.182806
#> epsilon (gene 1): -0.611677
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 3280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 562.356
#> delta (gene 1): 0.740415
#> phi (cell 1): 0.639178
#> s (cell 1): 1.5324
#> nu (cell 1): 0.440003
#> theta (batch 1): 0.415101
#> betas:    0.6780  -0.0527   0.2563   1.1200   0.4773   0.8416   0.5365   0.2696  -0.3344   0.1443  -0.1941  -0.4782
#> 
#> sigma: 0.182062
#> lambda (gene 1): 0.782203
#> epsilon (gene 1): -0.320629
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 3300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 593.252
#> delta (gene 1): 0.653455
#> phi (cell 1): 0.734051
#> s (cell 1): 0.849671
#> nu (cell 1): 0.424306
#> theta (batch 1): 0.301216
#> betas:    0.6550  -0.0465   0.0509   0.5534   0.7089   0.9733   0.2237   0.5839  -0.1384  -0.3803  -0.3180   0.0348
#> 
#> sigma: 0.176378
#> lambda (gene 1): 1.29649
#> epsilon (gene 1): -0.244199
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 3320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 809.359
#> delta (gene 1): 0.68642
#> phi (cell 1): 0.743572
#> s (cell 1): 0.213852
#> nu (cell 1): 0.407851
#> theta (batch 1): 0.548812
#> betas:    0.9935  -0.1072   0.3891   0.5296   0.2680   1.2885  -0.1036   0.5232   0.0603  -0.2969  -0.6775   0.2488
#> 
#> sigma: 0.167075
#> lambda (gene 1): 0.641204
#> epsilon (gene 1): -0.138813
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 3340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 642.736
#> delta (gene 1): 0.72949
#> phi (cell 1): 0.889665
#> s (cell 1): 0.164145
#> nu (cell 1): 0.408519
#> theta (batch 1): 0.652739
#> betas:    0.5595  -0.0464   0.2758   0.0347   1.0885   0.9751   0.3800   0.4770  -0.1165  -0.1965  -0.2961  -0.2139
#> 
#> sigma: 0.120607
#> lambda (gene 1): 1.44987
#> epsilon (gene 1): -0.171941
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 3360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 678.967
#> delta (gene 1): 0.586859
#> phi (cell 1): 0.959946
#> s (cell 1): 1.31411
#> nu (cell 1): 0.427641
#> theta (batch 1): 0.222303
#> betas:    0.9114  -0.1137   0.5204   0.5079   0.4487   0.7664   0.5245   0.1781  -0.1804   0.3132  -0.6644  -0.1172
#> 
#> sigma: 0.131213
#> lambda (gene 1): 1.33586
#> epsilon (gene 1): -0.348128
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -0.7
#> -----------------------------------------------------
#> Iteration 3380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 686.717
#> delta (gene 1): 0.655626
#> phi (cell 1): 0.753571
#> s (cell 1): 0.781122
#> nu (cell 1): 0.416718
#> theta (batch 1): 0.15733
#> betas:    0.7977  -0.0038   0.9149   0.8680   0.5326   0.7785   0.0753   0.3341   0.1319  -0.5330  -0.4068  -0.7132
#> 
#> sigma: 0.159396
#> lambda (gene 1): 0.513923
#> epsilon (gene 1): -0.232162
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.21822
#> LStheta (batch 1): -0.7
#> -----------------------------------------------------
#> Iteration 3400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 843.041
#> delta (gene 1): 0.571117
#> phi (cell 1): 0.772869
#> s (cell 1): 0.288466
#> nu (cell 1): 0.421146
#> theta (batch 1): 0.100431
#> betas:    0.6356  -0.1257   1.3913  -0.2602   1.0848   0.6992   0.5382   0.5927   0.2938   0.0036  -0.4502   0.3435
#> 
#> sigma: 0.231536
#> lambda (gene 1): 1.0068
#> epsilon (gene 1): -0.273707
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 6.6
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 3420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 692.958
#> delta (gene 1): 0.884255
#> phi (cell 1): 0.725022
#> s (cell 1): 1.84062
#> nu (cell 1): 0.447395
#> theta (batch 1): 0.213482
#> betas:    1.6427  -0.2911   0.5111   0.6216  -0.0456   0.2712   0.8346   0.4201  -0.0300  -0.1678   0.1217   0.0577
#> 
#> sigma: 0.171829
#> lambda (gene 1): 0.885555
#> epsilon (gene 1): 0.183172
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 6.6
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 3440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 432.472
#> delta (gene 1): 1.06001
#> phi (cell 1): 0.638726
#> s (cell 1): 1.13943
#> nu (cell 1): 0.416143
#> theta (batch 1): 0.269635
#> betas:    0.6287  -0.0788  -0.0213   0.8020   0.6124   0.8011   0.7953   0.1961   0.1005   0.0969  -0.3482  -0.3183
#> 
#> sigma: 0.201207
#> lambda (gene 1): 0.27778
#> epsilon (gene 1): 0.304116
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.2
#> LSphi: 6.6
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 3460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 762.885
#> delta (gene 1): 0.725564
#> phi (cell 1): 0.66726
#> s (cell 1): 0.289312
#> nu (cell 1): 0.422147
#> theta (batch 1): 0.0510394
#> betas:    1.0587  -0.1393   0.8390   0.7588   0.4792   0.3125   0.6429  -0.1811   0.1670   0.0370  -0.6683   0.1182
#> 
#> sigma: 0.15618
#> lambda (gene 1): 1.73125
#> epsilon (gene 1): -0.189922
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.7
#> -----------------------------------------------------
#> Iteration 3480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 750.09
#> delta (gene 1): 0.684279
#> phi (cell 1): 0.74218
#> s (cell 1): 0.414138
#> nu (cell 1): 0.423492
#> theta (batch 1): 0.0931346
#> betas:    1.1079  -0.1814  -0.2928   1.0907   0.3523   0.8252   0.1167   0.4613  -0.0466   0.0700  -0.2514  -0.1790
#> 
#> sigma: 0.173738
#> lambda (gene 1): 1.81047
#> epsilon (gene 1): -0.265202
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.3
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.7
#> -----------------------------------------------------
#> Iteration 3500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 545.071
#> delta (gene 1): 0.73437
#> phi (cell 1): 0.743547
#> s (cell 1): 0.393172
#> nu (cell 1): 0.42879
#> theta (batch 1): 0.0492356
#> betas:    0.9618  -0.2455   0.1632   0.8785   0.9958   0.0709   0.8706   0.8135   0.3202  -0.6930   0.8541   0.1060
#> 
#> sigma: 0.174045
#> lambda (gene 1): 2.44892
#> epsilon (gene 1): -0.0881338
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -0.6
#> -----------------------------------------------------
#> Iteration 3520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 715.227
#> delta (gene 1): 0.639768
#> phi (cell 1): 0.747678
#> s (cell 1): 0.421992
#> nu (cell 1): 0.418409
#> theta (batch 1): 0.0896091
#> betas:    0.7588  -0.1063   0.8326   0.4837   0.8250   0.9566   0.1811   0.2065   0.5827  -0.3108  -0.2209   0.0653
#> 
#> sigma: 0.141096
#> lambda (gene 1): 1.63084
#> epsilon (gene 1): -0.27734
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -0.6
#> -----------------------------------------------------
#> Iteration 3540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 614.575
#> delta (gene 1): 0.59714
#> phi (cell 1): 0.601902
#> s (cell 1): 0.363988
#> nu (cell 1): 0.413567
#> theta (batch 1): 0.224017
#> betas:    1.1615  -0.0897   0.3053   0.6530   1.0035   0.7498  -0.2583   0.4962   0.0022  -0.3362   0.0776  -1.2086
#> 
#> sigma: 0.165648
#> lambda (gene 1): 0.541289
#> epsilon (gene 1): -0.234029
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -0.6
#> -----------------------------------------------------
#> Iteration 3560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 567.928
#> delta (gene 1): 0.505536
#> phi (cell 1): 0.651267
#> s (cell 1): 0.558427
#> nu (cell 1): 0.424365
#> theta (batch 1): 0.343042
#> betas:    1.4332  -0.2503   0.5540   0.9790   0.7806  -0.0004   0.8300  -0.0312   0.7105  -0.7183  -0.0476   0.2563
#> 
#> sigma: 0.171742
#> lambda (gene 1): 1.54713
#> epsilon (gene 1): -0.338784
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.7
#> -----------------------------------------------------
#> Iteration 3580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 674.529
#> delta (gene 1): 0.695806
#> phi (cell 1): 0.642553
#> s (cell 1): 1.4725
#> nu (cell 1): 0.424136
#> theta (batch 1): 0.287515
#> betas:    0.7758  -0.0200   0.7467   0.9040   0.7353   0.3433   0.5087  -0.0579   0.2174  -0.2109  -0.8325  -0.1630
#> 
#> sigma: 0.169551
#> lambda (gene 1): 0.316228
#> epsilon (gene 1): -0.324346
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.7
#> -----------------------------------------------------
#> Iteration 3600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 932.876
#> delta (gene 1): 0.887832
#> phi (cell 1): 0.710984
#> s (cell 1): 0.566997
#> nu (cell 1): 0.43405
#> theta (batch 1): 0.206927
#> betas:    1.5418  -0.2223   0.4545   0.2648   0.9238   0.4255   0.2952   0.1465   0.4671  -0.2314  -0.2355  -0.1037
#> 
#> sigma: 0.21182
#> lambda (gene 1): 0.930516
#> epsilon (gene 1): 0.199451
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 3620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 778.078
#> delta (gene 1): 0.658054
#> phi (cell 1): 0.757371
#> s (cell 1): 0.434485
#> nu (cell 1): 0.431443
#> theta (batch 1): 0.115004
#> betas:    0.4403  -0.0371   0.0745  -0.0353   0.9377   0.9563   0.6784   0.3916  -0.0009   0.3696  -0.6415  -0.3118
#> 
#> sigma: 0.16937
#> lambda (gene 1): 1.56066
#> epsilon (gene 1): -0.0505105
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 3640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 587.035
#> delta (gene 1): 0.555728
#> phi (cell 1): 0.776324
#> s (cell 1): 0.32231
#> nu (cell 1): 0.42418
#> theta (batch 1): 0.128539
#> betas:    1.5004  -0.2883  -0.1708   1.1621   0.2954   1.3606  -0.0575   0.5207   0.3115   0.1948  -0.4302   0.4778
#> 
#> sigma: 0.168612
#> lambda (gene 1): 1.72167
#> epsilon (gene 1): -0.264167
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -0.8
#> -----------------------------------------------------
#> Iteration 3660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 672.958
#> delta (gene 1): 0.710404
#> phi (cell 1): 0.748135
#> s (cell 1): 0.568249
#> nu (cell 1): 0.429466
#> theta (batch 1): 0.176057
#> betas:    0.6442  -0.0636   0.1688   0.7056   0.3379   0.8016   0.3268   0.3466   0.0708  -0.1630  -0.0040  -0.5997
#> 
#> sigma: 0.149022
#> lambda (gene 1): 0.685106
#> epsilon (gene 1): -0.34888
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 3680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 769.863
#> delta (gene 1): 0.775041
#> phi (cell 1): 0.757617
#> s (cell 1): 0.697517
#> nu (cell 1): 0.423185
#> theta (batch 1): 0.05635
#> betas:    0.9390  -0.0507  -0.4434   1.0813   1.0605   0.1217   0.4094   0.5171  -0.5849  -0.1345  -0.0831  -0.8599
#> 
#> sigma: 0.228057
#> lambda (gene 1): 0.614836
#> epsilon (gene 1): -0.150167
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 3700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 483.974
#> delta (gene 1): 0.956729
#> phi (cell 1): 0.679997
#> s (cell 1): 0.474161
#> nu (cell 1): 0.419476
#> theta (batch 1): 0.118678
#> betas:    0.7386  -0.1638   0.5793  -0.1724   1.0971   0.9252   0.6347   0.1671   0.1128   0.2057  -0.1233   0.0975
#> 
#> sigma: 0.210681
#> lambda (gene 1): 0.353458
#> epsilon (gene 1): 0.000894185
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 752.457
#> delta (gene 1): 0.719839
#> phi (cell 1): 0.803728
#> s (cell 1): 0.406185
#> nu (cell 1): 0.40993
#> theta (batch 1): 0.116192
#> betas:    0.9836  -0.0480   0.8044  -0.0251   1.1700   0.3023   0.4344  -0.3416   0.2485  -0.4683  -0.5728  -0.1541
#> 
#> sigma: 0.172934
#> lambda (gene 1): 1.3801
#> epsilon (gene 1): -0.161047
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 676.021
#> delta (gene 1): 0.655025
#> phi (cell 1): 0.735517
#> s (cell 1): 0.552083
#> nu (cell 1): 0.415364
#> theta (batch 1): 0.202655
#> betas:    0.3362   0.0120   0.5488   1.1026   0.6518   0.5111   0.6023   0.0606   0.1046  -0.2859  -0.4951  -0.2354
#> 
#> sigma: 0.167631
#> lambda (gene 1): 2.58947
#> epsilon (gene 1): -0.20977
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 3760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 746.504
#> delta (gene 1): 1.10915
#> phi (cell 1): 0.61936
#> s (cell 1): 0.380326
#> nu (cell 1): 0.426977
#> theta (batch 1): 0.174002
#> betas:    1.0921  -0.1159  -0.0275   1.0491   0.7821   0.2025   0.1015   0.4290  -0.1019  -0.2196  -0.2596  -0.3327
#> 
#> sigma: 0.170242
#> lambda (gene 1): 0.440287
#> epsilon (gene 1): 0.146448
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 3780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 546.35
#> delta (gene 1): 0.785154
#> phi (cell 1): 0.605609
#> s (cell 1): 0.177884
#> nu (cell 1): 0.439653
#> theta (batch 1): 0.376474
#> betas:    1.0240  -0.0992   0.4929   0.8455   0.5900   0.5045   0.0165   0.0672   0.3675  -0.6125   0.1247  -0.5887
#> 
#> sigma: 0.137441
#> lambda (gene 1): 1.32649
#> epsilon (gene 1): -0.290667
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 3800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 712.622
#> delta (gene 1): 0.789572
#> phi (cell 1): 0.64017
#> s (cell 1): 0.374489
#> nu (cell 1): 0.417783
#> theta (batch 1): 0.320065
#> betas:    0.9944  -0.1031   0.7278   0.1915   0.7804   0.7105   0.3649  -0.0700   0.5517  -0.4067  -0.1080  -0.9075
#> 
#> sigma: 0.183549
#> lambda (gene 1): 1.54531
#> epsilon (gene 1): -0.00267522
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 3820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 696.329
#> delta (gene 1): 1.05304
#> phi (cell 1): 0.646971
#> s (cell 1): 2.74042
#> nu (cell 1): 0.411766
#> theta (batch 1): 0.37311
#> betas:    0.7902  -0.1060   0.6195   0.9324   0.5908   0.7089   0.5570   0.0025   0.0190  -0.2163  -0.0434  -0.2248
#> 
#> sigma: 0.163095
#> lambda (gene 1): 1.13496
#> epsilon (gene 1): 0.271909
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 3840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 502.114
#> delta (gene 1): 0.84491
#> phi (cell 1): 0.631803
#> s (cell 1): 0.845461
#> nu (cell 1): 0.414417
#> theta (batch 1): 0.417122
#> betas:    1.2714  -0.1448   0.0491   1.2900   0.1801   0.4193   0.3136   0.5260  -0.2635  -0.4106  -0.0086  -0.2636
#> 
#> sigma: 0.130283
#> lambda (gene 1): 0.943261
#> epsilon (gene 1): 0.0762994
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.4
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 3860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 668.559
#> delta (gene 1): 0.746354
#> phi (cell 1): 0.67938
#> s (cell 1): 0.47927
#> nu (cell 1): 0.419428
#> theta (batch 1): 0.516508
#> betas:    0.7917  -0.0614   0.2278  -0.1118   1.4120   0.6087   0.3888   0.2787   0.0473  -0.3292  -0.4135  -0.2767
#> 
#> sigma: 0.142173
#> lambda (gene 1): 0.798641
#> epsilon (gene 1): 0.101268
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 3880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 682.592
#> delta (gene 1): 0.587927
#> phi (cell 1): 0.738066
#> s (cell 1): 1.50929
#> nu (cell 1): 0.434147
#> theta (batch 1): 0.778464
#> betas:    0.8110  -0.0887   0.6973  -0.3237   0.8593   0.8770   0.6170   0.1767   0.2895  -0.3600  -0.1772  -0.1154
#> 
#> sigma: 0.118954
#> lambda (gene 1): 0.441087
#> epsilon (gene 1): -0.340955
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 3900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1006.59
#> delta (gene 1): 0.675781
#> phi (cell 1): 0.733081
#> s (cell 1): 0.71602
#> nu (cell 1): 0.425361
#> theta (batch 1): 0.224515
#> betas:    1.1541  -0.2588  -0.2169   1.0822   0.1838   1.1628   0.3905   0.1643   0.5504   0.0635  -0.0569   0.2081
#> 
#> sigma: 0.22217
#> lambda (gene 1): 1.33249
#> epsilon (gene 1): -0.111268
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 3920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 992.143
#> delta (gene 1): 0.86563
#> phi (cell 1): 0.911064
#> s (cell 1): 0.420771
#> nu (cell 1): 0.415438
#> theta (batch 1): 0.0658736
#> betas:    0.4519   0.0604   0.6022   0.6894   1.0077   0.1295   0.4851   0.5039  -0.4343  -0.7431   0.1348  -1.1062
#> 
#> sigma: 0.16277
#> lambda (gene 1): 1.51571
#> epsilon (gene 1): 0.199129
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 3940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 658.167
#> delta (gene 1): 0.591427
#> phi (cell 1): 0.872999
#> s (cell 1): 0.3937
#> nu (cell 1): 0.428531
#> theta (batch 1): 0.0899172
#> betas:    1.2099  -0.1861   0.4689   0.2570   0.6241   0.7602   0.2915   0.1393   0.2239  -0.2505  -0.1509  -0.1606
#> 
#> sigma: 0.140481
#> lambda (gene 1): 0.80835
#> epsilon (gene 1): -0.112298
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.4
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 3960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 628.593
#> delta (gene 1): 0.659433
#> phi (cell 1): 0.783649
#> s (cell 1): 0.498098
#> nu (cell 1): 0.409282
#> theta (batch 1): 0.120497
#> betas:    0.1076  -0.0047   0.7979   1.2800   0.8457   0.9290   0.5539   0.4052   0.2049  -0.1983  -0.1155  -0.4455
#> 
#> sigma: 0.15763
#> lambda (gene 1): 0.918532
#> epsilon (gene 1): -0.318273
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 3980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 633.37
#> delta (gene 1): 0.903508
#> phi (cell 1): 0.83491
#> s (cell 1): 0.351845
#> nu (cell 1): 0.428847
#> theta (batch 1): 0.124189
#> betas:    0.8769   0.0283  -0.1320   0.3243   0.5810   0.7288   0.3620  -0.0591  -0.0953  -0.0972  -0.7406  -1.2707
#> 
#> sigma: 0.14412
#> lambda (gene 1): 0.485534
#> epsilon (gene 1): -0.0643368
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 4000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 725.343
#> delta (gene 1): 0.718594
#> phi (cell 1): 0.866288
#> s (cell 1): 0.457554
#> nu (cell 1): 0.41995
#> theta (batch 1): 0.0504906
#> betas:    0.5492  -0.0862   0.3453   0.7951   0.3343   0.9619   0.9091   0.0909   0.1489   0.1913  -0.6534  -0.0526
#> 
#> sigma: 0.17428
#> lambda (gene 1): 1.24586
#> epsilon (gene 1): -0.00780468
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 4020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 708.47
#> delta (gene 1): 0.85376
#> phi (cell 1): 0.890365
#> s (cell 1): 0.42388
#> nu (cell 1): 0.425852
#> theta (batch 1): 0.0620792
#> betas:    1.1151  -0.1582   0.3940   0.8461  -0.2505   0.8526   0.7578  -0.1555   0.1304   0.2977  -0.2233  -0.5211
#> 
#> sigma: 0.182192
#> lambda (gene 1): 0.699788
#> epsilon (gene 1): -0.16939
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 4040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 649.423
#> delta (gene 1): 1.00152
#> phi (cell 1): 0.694781
#> s (cell 1): 0.762385
#> nu (cell 1): 0.418426
#> theta (batch 1): 0.0599253
#> betas:    0.6904  -0.0258   0.6656  -0.0749   0.6624   1.1054   0.7676  -0.0809   0.0336  -0.2848   0.0572  -1.0540
#> 
#> sigma: 0.192872
#> lambda (gene 1): 0.619508
#> epsilon (gene 1): 0.104431
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.6
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 4060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 818.476
#> delta (gene 1): 0.808487
#> phi (cell 1): 0.610426
#> s (cell 1): 0.409183
#> nu (cell 1): 0.423121
#> theta (batch 1): 0.00598497
#> betas:    0.2256  -0.0296   0.4904   0.6360   1.1627   0.5949   0.6801   0.1839   0.1041  -0.1206  -0.0687  -0.4838
#> 
#> sigma: 0.152765
#> lambda (gene 1): 1.05919
#> epsilon (gene 1): -0.0514676
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 4080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 718.356
#> delta (gene 1): 0.658633
#> phi (cell 1): 0.628565
#> s (cell 1): 0.435509
#> nu (cell 1): 0.437416
#> theta (batch 1): 0.00952505
#> betas:    1.3339  -0.1810   0.6047  -0.0661   0.9167   0.6033  -0.1956   0.8416  -0.2820   0.0128  -0.3210  -0.1444
#> 
#> sigma: 0.157722
#> lambda (gene 1): 0.920211
#> epsilon (gene 1): -0.277783
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.5
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 4100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 797.701
#> delta (gene 1): 0.686465
#> phi (cell 1): 0.647772
#> s (cell 1): 0.401323
#> nu (cell 1): 0.41944
#> theta (batch 1): 0.0105714
#> betas:    0.5427  -0.0653   0.5938   0.5101   1.4426   0.2552   0.6440   0.3476  -0.0662   0.0391  -0.2754  -0.5128
#> 
#> sigma: 0.166923
#> lambda (gene 1): 0.518116
#> epsilon (gene 1): -0.102684
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 725.435
#> delta (gene 1): 0.552704
#> phi (cell 1): 0.702438
#> s (cell 1): 0.383938
#> nu (cell 1): 0.429017
#> theta (batch 1): 0.0104424
#> betas:    0.7137  -0.1438   0.7588   0.6310   1.0338   0.4857   0.4113   0.1082   0.5292  -0.4206   0.1166  -0.0694
#> 
#> sigma: 0.127083
#> lambda (gene 1): 0.489159
#> epsilon (gene 1): -0.313968
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 795.054
#> delta (gene 1): 0.861534
#> phi (cell 1): 0.662719
#> s (cell 1): 0.441901
#> nu (cell 1): 0.431791
#> theta (batch 1): 0.0103163
#> betas:    0.9111  -0.0145   0.1661   0.6718   0.9372   0.3806   0.0605  -0.2429   0.3004  -0.3413  -0.8068  -0.4686
#> 
#> sigma: 0.162302
#> lambda (gene 1): 1.75461
#> epsilon (gene 1): 0.126586
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 606.744
#> delta (gene 1): 0.786307
#> phi (cell 1): 0.707368
#> s (cell 1): 0.423858
#> nu (cell 1): 0.413265
#> theta (batch 1): 0.00873159
#> betas:    0.5728  -0.0032  -0.0716   0.5779   1.0168   0.3335   0.4926   0.5052  -0.1314  -0.2256  -0.4609  -0.1687
#> 
#> sigma: 0.177704
#> lambda (gene 1): 0.694398
#> epsilon (gene 1): -0.11314
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 4180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 621.071
#> delta (gene 1): 0.600092
#> phi (cell 1): 0.614403
#> s (cell 1): 0.395842
#> nu (cell 1): 0.433565
#> theta (batch 1): 0.00482343
#> betas:    1.3184  -0.1547   0.0371  -0.0134   0.4778   1.0787   0.3540   0.0478   0.0100   0.1379  -0.6483  -0.2751
#> 
#> sigma: 0.154774
#> lambda (gene 1): 0.230548
#> epsilon (gene 1): -0.274729
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 4200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 766.444
#> delta (gene 1): 0.814043
#> phi (cell 1): 0.644869
#> s (cell 1): 0.484705
#> nu (cell 1): 0.423446
#> theta (batch 1): 0.0362836
#> betas:    1.0785  -0.1144   0.4901   0.5833   0.4765   0.7652   0.2510   0.3908  -0.3723  -0.0251  -0.4006  -0.3170
#> 
#> sigma: 0.165615
#> lambda (gene 1): 2.40004
#> epsilon (gene 1): 0.101285
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.8
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 4220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 734.444
#> delta (gene 1): 0.518726
#> phi (cell 1): 0.674777
#> s (cell 1): 0.351191
#> nu (cell 1): 0.421487
#> theta (batch 1): 0.102918
#> betas:    0.5960  -0.0304   0.2645   1.0841   0.8302   0.3672   0.8862  -0.0066  -0.0312   0.0377  -0.5944  -0.0973
#> 
#> sigma: 0.164241
#> lambda (gene 1): 0.958922
#> epsilon (gene 1): -0.538597
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.8
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 4240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 881.222
#> delta (gene 1): 0.736198
#> phi (cell 1): 0.841469
#> s (cell 1): 0.580877
#> nu (cell 1): 0.422617
#> theta (batch 1): 0.307245
#> betas:    0.5928  -0.2524   0.5360   0.4038   0.8109   0.6923   0.6705   1.0144   0.1768   0.1102   0.1082   0.5873
#> 
#> sigma: 0.203737
#> lambda (gene 1): 2.579
#> epsilon (gene 1): 0.0919383
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.8
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 4260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 720.137
#> delta (gene 1): 0.691873
#> phi (cell 1): 0.73216
#> s (cell 1): 0.862848
#> nu (cell 1): 0.440281
#> theta (batch 1): 0.101237
#> betas:    1.2020  -0.1560  -0.2494   1.1273   0.5234   0.7387   0.2849   0.1544   0.1576  -0.0705  -0.5238   0.1923
#> 
#> sigma: 0.147985
#> lambda (gene 1): 0.383634
#> epsilon (gene 1): -0.287997
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.5
#> -----------------------------------------------------
#> Iteration 4280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 788.645
#> delta (gene 1): 0.802758
#> phi (cell 1): 0.927962
#> s (cell 1): 0.441492
#> nu (cell 1): 0.405519
#> theta (batch 1): 0.101809
#> betas:    0.2620  -0.0006   0.8875   0.4092   1.4899   0.6828   0.1639   0.4443  -0.1170   0.2772  -1.0930   0.3865
#> 
#> sigma: 0.136865
#> lambda (gene 1): 0.27183
#> epsilon (gene 1): -0.0600197
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.5
#> -----------------------------------------------------
#> Iteration 4300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 843.338
#> delta (gene 1): 0.907026
#> phi (cell 1): 1.08274
#> s (cell 1): 0.967164
#> nu (cell 1): 0.432876
#> theta (batch 1): 0.10806
#> betas:    0.6138   0.1071   0.4445   0.6257   0.6312   0.5812   0.3722  -0.1395  -0.4401   0.0078  -0.8604  -1.1206
#> 
#> sigma: 0.17613
#> lambda (gene 1): 0.981663
#> epsilon (gene 1): 0.0775536
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.8
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.6
#> -----------------------------------------------------
#> Iteration 4320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 835.051
#> delta (gene 1): 0.591295
#> phi (cell 1): 0.999356
#> s (cell 1): 1.21051
#> nu (cell 1): 0.418027
#> theta (batch 1): 0.445399
#> betas:    0.6418  -0.0440   0.1944   0.7507   0.6482   1.0959   0.1374   0.2822   0.0553   0.1563  -0.7277  -0.1395
#> 
#> sigma: 0.167466
#> lambda (gene 1): 0.839534
#> epsilon (gene 1): -0.364121
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.8
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.6
#> -----------------------------------------------------
#> Iteration 4340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 648.071
#> delta (gene 1): 0.812483
#> phi (cell 1): 0.944854
#> s (cell 1): 0.428239
#> nu (cell 1): 0.419793
#> theta (batch 1): 0.319393
#> betas:    0.9408  -0.0311   0.3447   1.2043   0.9437   0.6448  -0.0166   0.3663   0.0023  -0.4241  -0.4348  -0.5979
#> 
#> sigma: 0.151721
#> lambda (gene 1): 2.12711
#> epsilon (gene 1): 0.0626462
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.8
#> LSphi: 7
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.6
#> -----------------------------------------------------
#> Iteration 4360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 713.91
#> delta (gene 1): 0.845871
#> phi (cell 1): 1.00055
#> s (cell 1): 1.90092
#> nu (cell 1): 0.410774
#> theta (batch 1): 0.381544
#> betas:    0.8991  -0.1238   0.6638   0.2120   0.8262   0.8648   0.2932   0.1646   0.0805  -0.0364  -0.2319  -0.0574
#> 
#> sigma: 0.176699
#> lambda (gene 1): 1.31248
#> epsilon (gene 1): -0.027137
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.9
#> LSphi: 6.9
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -1.5
#> -----------------------------------------------------
#> Iteration 4380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 885.454
#> delta (gene 1): 0.635609
#> phi (cell 1): 1.12462
#> s (cell 1): 0.458379
#> nu (cell 1): 0.425052
#> theta (batch 1): 0.120652
#> betas:    0.8008  -0.1626   1.0055   0.6499   0.8021   0.9148   0.1812  -0.0536   1.0236  -0.4600  -0.1654   0.1720
#> 
#> sigma: 0.204328
#> lambda (gene 1): 2.60342
#> epsilon (gene 1): -0.144235
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.9
#> LSphi: 6.9
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -1.5
#> -----------------------------------------------------
#> Iteration 4400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 667.847
#> delta (gene 1): 1.1007
#> phi (cell 1): 0.991405
#> s (cell 1): 0.290474
#> nu (cell 1): 0.426992
#> theta (batch 1): 0.200875
#> betas:    0.9289  -0.0888   0.3026   0.6288   0.8596   0.8225  -0.0617   0.3162   0.3202  -0.5564   0.1056  -0.7439
#> 
#> sigma: 0.183159
#> lambda (gene 1): 2.43853
#> epsilon (gene 1): 0.203564
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.8
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 4420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 986.096
#> delta (gene 1): 0.747513
#> phi (cell 1): 0.965092
#> s (cell 1): 0.563242
#> nu (cell 1): 0.417452
#> theta (batch 1): 0.227347
#> betas:    0.5743  -0.1165   0.5748   0.8705   1.1568   0.8849   0.1542   0.6443  -0.0494   0.0377   0.2070  -0.3396
#> 
#> sigma: 0.162335
#> lambda (gene 1): 1.68805
#> epsilon (gene 1): -0.296337
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.8
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 4440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 867.685
#> delta (gene 1): 0.892705
#> phi (cell 1): 0.819762
#> s (cell 1): 1.35869
#> nu (cell 1): 0.424112
#> theta (batch 1): 0.138961
#> betas:    0.9934  -0.1559   0.4531   0.3283   0.6956   0.7790   0.5735  -0.0230   0.3633  -0.0194  -0.2595  -0.1121
#> 
#> sigma: 0.147165
#> lambda (gene 1): 0.724212
#> epsilon (gene 1): -0.061013
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.8
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.4
#> -----------------------------------------------------
#> Iteration 4460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 713.328
#> delta (gene 1): 0.825272
#> phi (cell 1): 0.749232
#> s (cell 1): 0.459173
#> nu (cell 1): 0.407278
#> theta (batch 1): 0.238782
#> betas:    0.8521  -0.1144   0.8901   0.3890   0.6077   0.8461   0.3152   0.1855   0.0326  -0.1244  -0.3766  -0.0369
#> 
#> sigma: 0.159288
#> lambda (gene 1): 0.550342
#> epsilon (gene 1): -0.09343
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.9
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 4480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 666.698
#> delta (gene 1): 0.433832
#> phi (cell 1): 0.74255
#> s (cell 1): 0.217154
#> nu (cell 1): 0.427908
#> theta (batch 1): 0.191091
#> betas:    0.2972   0.0494   0.7629   0.1030   0.9669   1.1475  -0.5880   0.2350   0.5117  -0.8963  -0.6471  -0.0693
#> 
#> sigma: 0.146314
#> lambda (gene 1): 1.34369
#> epsilon (gene 1): -0.659977
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.9
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.3
#> -----------------------------------------------------
#> Iteration 4500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 777.887
#> delta (gene 1): 0.513754
#> phi (cell 1): 0.784089
#> s (cell 1): 0.29146
#> nu (cell 1): 0.422768
#> theta (batch 1): 0.471864
#> betas:    0.7950  -0.0514   0.4216   0.9798   0.4747   0.6165   0.4348   0.4109  -0.4716  -0.3634   0.1484  -0.8591
#> 
#> sigma: 0.152746
#> lambda (gene 1): 0.247721
#> epsilon (gene 1): -0.411828
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -2
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 755.728
#> delta (gene 1): 0.735696
#> phi (cell 1): 0.711444
#> s (cell 1): 0.360947
#> nu (cell 1): 0.434953
#> theta (batch 1): 0.286442
#> betas:    0.8178  -0.0772   0.5825   0.4449   0.7097   0.8543  -0.0657   0.6812  -0.1813  -0.1749  -0.2432  -0.4488
#> 
#> sigma: 0.177264
#> lambda (gene 1): 0.359176
#> epsilon (gene 1): -0.125606
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -2
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 656.588
#> delta (gene 1): 0.527193
#> phi (cell 1): 0.715976
#> s (cell 1): 0.760389
#> nu (cell 1): 0.418675
#> theta (batch 1): 0.349981
#> betas:    0.2432   0.1086   0.2061   0.6749   0.3242   1.3124   0.1462   0.3667  -0.5134  -0.0751  -0.7870  -0.6703
#> 
#> sigma: 0.172051
#> lambda (gene 1): 0.657169
#> epsilon (gene 1): -0.428936
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -2
#> LSphi: 6.8
#> LSnu (cell 1): -6.11822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 697.015
#> delta (gene 1): 0.60156
#> phi (cell 1): 0.869968
#> s (cell 1): 0.240639
#> nu (cell 1): 0.39725
#> theta (batch 1): 0.211811
#> betas:    0.7286  -0.0972   0.2085   0.6886   0.9105   1.0645  -0.1663   0.7771  -0.2651   0.1506  -0.4082  -0.0979
#> 
#> sigma: 0.18474
#> lambda (gene 1): 0.941324
#> epsilon (gene 1): -0.291351
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.9
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 4580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 675.84
#> delta (gene 1): 0.652204
#> phi (cell 1): 0.837907
#> s (cell 1): 0.395996
#> nu (cell 1): 0.429429
#> theta (batch 1): 0.256247
#> betas:    1.0179  -0.1101  -0.3191   0.4674   0.6546   0.7800   0.6013   0.0541  -0.0650  -0.1004  -0.3342  -0.2976
#> 
#> sigma: 0.136481
#> lambda (gene 1): 0.771807
#> epsilon (gene 1): -0.180262
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.9
#> LSphi: 6.7
#> LSnu (cell 1): -6.01822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 4600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1223.49
#> delta (gene 1): 0.932906
#> phi (cell 1): 0.723173
#> s (cell 1): 0.327943
#> nu (cell 1): 0.41676
#> theta (batch 1): 0.0593314
#> betas:    0.0390  -0.0348   0.6552   0.7605   0.7265   1.2150   0.3623   0.4714   0.4480  -0.0917  -0.3450  -0.0689
#> 
#> sigma: 0.161482
#> lambda (gene 1): 0.977735
#> epsilon (gene 1): 0.289268
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -2
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 712.277
#> delta (gene 1): 0.534784
#> phi (cell 1): 0.700487
#> s (cell 1): 0.362054
#> nu (cell 1): 0.424742
#> theta (batch 1): 0.0117024
#> betas:    0.7306  -0.0588   0.0785   0.2898   0.7195   0.7970   0.8849  -0.1146  -0.0696  -0.0293  -0.4648  -0.3430
#> 
#> sigma: 0.14928
#> lambda (gene 1): 1.32215
#> epsilon (gene 1): -0.317651
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -2
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 837.404
#> delta (gene 1): 0.973433
#> phi (cell 1): 0.706602
#> s (cell 1): 0.534666
#> nu (cell 1): 0.419246
#> theta (batch 1): 0.043631
#> betas:    0.4289  -0.0594   0.8448   0.6595   0.6628   0.9748   0.7800   0.2655  -0.1545   0.3108  -0.4974  -0.2345
#> 
#> sigma: 0.198599
#> lambda (gene 1): 0.0819747
#> epsilon (gene 1): 0.250657
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -2
#> LSphi: 6.8
#> LSnu (cell 1): -5.91822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 915.945
#> delta (gene 1): 0.677754
#> phi (cell 1): 0.761253
#> s (cell 1): 0.449457
#> nu (cell 1): 0.410267
#> theta (batch 1): 0.0854119
#> betas:    0.9816  -0.0375   0.6344   0.5115   0.5296   0.6458  -0.0142  -0.0258   0.2599  -0.2784  -0.6471  -0.6046
#> 
#> sigma: 0.171874
#> lambda (gene 1): 2.36632
#> epsilon (gene 1): -0.171199
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.9
#> LSphi: 6.7
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 4680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 571.731
#> delta (gene 1): 0.656272
#> phi (cell 1): 0.677595
#> s (cell 1): 0.579992
#> nu (cell 1): 0.434144
#> theta (batch 1): 0.204597
#> betas:    0.2611   0.0486   1.2404   0.3277   0.5959   1.4528   0.0940   0.1334  -0.3108   0.2739  -0.7666  -0.4413
#> 
#> sigma: 0.223601
#> lambda (gene 1): 0.560088
#> epsilon (gene 1): -0.272886
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.9
#> LSphi: 6.7
#> LSnu (cell 1): -5.81822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 4700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 590.239
#> delta (gene 1): 0.82475
#> phi (cell 1): 0.650378
#> s (cell 1): 2.74826
#> nu (cell 1): 0.410483
#> theta (batch 1): 0.510654
#> betas:    0.4591  -0.0174   1.1323   0.1645   1.0308   0.5725   0.9284  -0.0095   0.5013  -0.0546  -0.7271  -0.1458
#> 
#> sigma: 0.14935
#> lambda (gene 1): 0.915392
#> epsilon (gene 1): -0.0659657
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.8
#> LSphi: 6.8
#> LSnu (cell 1): -5.71822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 650.594
#> delta (gene 1): 0.707496
#> phi (cell 1): 0.654921
#> s (cell 1): 0.296139
#> nu (cell 1): 0.430434
#> theta (batch 1): 0.302287
#> betas:    1.2004  -0.0205   0.3041   0.1458   0.2635   0.8512   0.1916   0.3000  -0.2982  -0.2204  -0.4313  -1.2225
#> 
#> sigma: 0.146716
#> lambda (gene 1): 1.40943
#> epsilon (gene 1): -0.226257
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.8
#> LSphi: 6.8
#> LSnu (cell 1): -5.71822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 798.183
#> delta (gene 1): 0.495234
#> phi (cell 1): 0.652334
#> s (cell 1): 0.253153
#> nu (cell 1): 0.427066
#> theta (batch 1): 0.260742
#> betas:    1.0882  -0.1683   0.3427   0.0153   0.8843   0.4419   0.8792  -0.1287   0.6370  -0.0833  -0.6793   0.2493
#> 
#> sigma: 0.164701
#> lambda (gene 1): 0.847866
#> epsilon (gene 1): -0.469041
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.2
#> LSdelta (gene 1): -1.8
#> LSphi: 6.8
#> LSnu (cell 1): -5.71822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 728.35
#> delta (gene 1): 0.570126
#> phi (cell 1): 0.656316
#> s (cell 1): 0.620918
#> nu (cell 1): 0.431748
#> theta (batch 1): 0.326574
#> betas:    0.6813   0.0145   0.4237   0.6846   0.7985   0.7862   0.0553  -0.0142   0.2007  -0.6782  -0.0496  -0.9855
#> 
#> sigma: 0.136819
#> lambda (gene 1): 0.343447
#> epsilon (gene 1): -0.305823
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -5.61822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 4780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 854.24
#> delta (gene 1): 1.00068
#> phi (cell 1): 0.710564
#> s (cell 1): 0.56233
#> nu (cell 1): 0.42773
#> theta (batch 1): 0.162663
#> betas:    1.2513  -0.1639   0.0895   0.3885   0.9206   0.4252   0.5418  -0.2172   0.0061   0.0881  -0.4422   0.1873
#> 
#> sigma: 0.143194
#> lambda (gene 1): 2.26192
#> epsilon (gene 1): 0.115178
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -5.61822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 4800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 894.419
#> delta (gene 1): 0.553485
#> phi (cell 1): 0.700942
#> s (cell 1): 0.649375
#> nu (cell 1): 0.421448
#> theta (batch 1): 0.198846
#> betas:    1.4592  -0.2047  -0.2561   0.6847   0.5837   0.4407   0.3478   0.3087  -0.0832   0.0400  -0.5852  -0.0070
#> 
#> sigma: 0.110219
#> lambda (gene 1): 1.96812
#> epsilon (gene 1): -0.11288
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.8
#> LSphi: 7
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 581.415
#> delta (gene 1): 0.604692
#> phi (cell 1): 0.795619
#> s (cell 1): 0.925776
#> nu (cell 1): 0.420347
#> theta (batch 1): 0.35011
#> betas:    0.6091  -0.0558   0.1355   0.6641   0.2499   1.0580   0.3912  -0.0776   0.2596   0.0515  -0.8004  -0.0505
#> 
#> sigma: 0.153618
#> lambda (gene 1): 0.673192
#> epsilon (gene 1): -0.0808115
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.8
#> LSphi: 7
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 667.876
#> delta (gene 1): 0.691328
#> phi (cell 1): 0.687668
#> s (cell 1): 0.693303
#> nu (cell 1): 0.424919
#> theta (batch 1): 0.275754
#> betas:    0.8909  -0.1075   0.6859   0.4127   0.7812   0.5179   0.4620   0.3894  -0.3109   0.3191  -0.3120  -0.5898
#> 
#> sigma: 0.221071
#> lambda (gene 1): 1.65842
#> epsilon (gene 1): -0.106012
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.8
#> LSphi: 7
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1.2
#> -----------------------------------------------------
#> Iteration 4860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 688.93
#> delta (gene 1): 0.636738
#> phi (cell 1): 0.74263
#> s (cell 1): 0.591042
#> nu (cell 1): 0.438799
#> theta (batch 1): 0.159457
#> betas:    0.4422  -0.0105   0.3881   0.4262   1.1091   0.5885   0.4183   0.3661  -0.2152   0.1632  -0.7707  -0.0357
#> 
#> sigma: 0.141174
#> lambda (gene 1): 1.5303
#> epsilon (gene 1): -0.274628
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -5.61822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 4880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 756.612
#> delta (gene 1): 0.859348
#> phi (cell 1): 0.794348
#> s (cell 1): 0.332692
#> nu (cell 1): 0.420942
#> theta (batch 1): 0.318417
#> betas:    0.1669   0.1748   0.5437   0.6333   0.6334   0.5611   0.1039   0.0358   0.2946  -0.9288  -0.7719  -0.9692
#> 
#> sigma: 0.182117
#> lambda (gene 1): 0.69513
#> epsilon (gene 1): 0.250528
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.3
#> LSdelta (gene 1): -1.7
#> LSphi: 6.9
#> LSnu (cell 1): -5.61822
#> LStheta (batch 1): -1.1
#> -----------------------------------------------------
#> Iteration 4900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 719.131
#> delta (gene 1): 0.644975
#> phi (cell 1): 0.740866
#> s (cell 1): 0.721502
#> nu (cell 1): 0.407595
#> theta (batch 1): 0.365064
#> betas:    0.6525  -0.0437   1.1484  -0.2528   0.7168   1.0793   0.4518  -0.0887   0.2484  -0.0396  -0.3461  -0.4889
#> 
#> sigma: 0.15525
#> lambda (gene 1): 0.499325
#> epsilon (gene 1): -0.380381
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 4920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 698.034
#> delta (gene 1): 0.689192
#> phi (cell 1): 0.735165
#> s (cell 1): 1.1996
#> nu (cell 1): 0.425567
#> theta (batch 1): 0.219502
#> betas:    1.1961  -0.1970   0.1055   0.2161   0.4170   0.9678   0.4955   0.3318  -0.0619  -0.0192  -0.1408   0.0861
#> 
#> sigma: 0.199174
#> lambda (gene 1): 1.24248
#> epsilon (gene 1): -0.1889
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 4940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 695.153
#> delta (gene 1): 0.501943
#> phi (cell 1): 0.756578
#> s (cell 1): 0.520823
#> nu (cell 1): 0.428388
#> theta (batch 1): 0.0673934
#> betas:    0.9413  -0.1465   0.1456  -0.0652   0.6545   0.9432   0.7064   0.2577   0.0405   0.2223  -0.6259   0.1990
#> 
#> sigma: 0.203456
#> lambda (gene 1): 0.862084
#> epsilon (gene 1): -0.474365
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 4960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 627.854
#> delta (gene 1): 0.691076
#> phi (cell 1): 0.808759
#> s (cell 1): 0.471172
#> nu (cell 1): 0.412812
#> theta (batch 1): 0.0756824
#> betas:    1.2032  -0.2335   0.5983   0.8364   0.4661   0.7747   0.2584   0.2623  -0.0118   0.1787  -0.0031  -0.1802
#> 
#> sigma: 0.141623
#> lambda (gene 1): 0.893593
#> epsilon (gene 1): -0.239563
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.5
#> LSphi: 7.1
#> LSnu (cell 1): -5.41822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> Iteration 4980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 700.303
#> delta (gene 1): 0.700406
#> phi (cell 1): 0.882286
#> s (cell 1): 0.298498
#> nu (cell 1): 0.418171
#> theta (batch 1): 0.369819
#> betas:    1.3417  -0.1986   0.4276   0.2908   0.8793   0.7164   0.2540  -0.0646   0.2196   0.0923  -0.4606  -0.1706
#> 
#> sigma: 0.230809
#> lambda (gene 1): 1.57042
#> epsilon (gene 1): -0.252079
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.5
#> LSdelta (gene 1): -1.5
#> LSphi: 7.1
#> LSnu (cell 1): -5.41822
#> LStheta (batch 1): -0.9
#> -----------------------------------------------------
#> End of Burn-in period.
#> -----------------------------------------------------
#> -----------------------------------------------------
#> Iteration 5000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 495.147
#> delta (gene 1): 0.841352
#> phi (cell 1): 1.03287
#> s (cell 1): 2.07827
#> nu (cell 1): 0.441795
#> theta (batch 1): 0.331239
#> betas:    0.6921  -0.0162   0.6480   1.2028   0.1475   0.5975   0.4784  -0.3898   0.6508  -0.9667   0.1207  -0.6750
#> 
#> sigma: 0.18332
#> lambda (gene 1): 0.426813
#> epsilon (gene 1): -0.047719
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 724.619
#> delta (gene 1): 0.769575
#> phi (cell 1): 1.02011
#> s (cell 1): 0.67217
#> nu (cell 1): 0.424205
#> theta (batch 1): 0.296361
#> betas:    1.7567  -0.3517   0.8302   0.2184   0.6250   0.5428   0.4094   0.2868   0.2377  -0.3415   0.3194   0.4111
#> 
#> sigma: 0.138534
#> lambda (gene 1): 2.82001
#> epsilon (gene 1): 0.0677743
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 566.327
#> delta (gene 1): 0.716932
#> phi (cell 1): 1.01167
#> s (cell 1): 0.412197
#> nu (cell 1): 0.418183
#> theta (batch 1): 0.637074
#> betas:    0.3898  -0.0803   1.2810   0.9595   0.5985   0.3811   0.6037   0.2899   0.1361   0.0763  -0.3572   0.1384
#> 
#> sigma: 0.149439
#> lambda (gene 1): 0.976539
#> epsilon (gene 1): -0.203803
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 612.056
#> delta (gene 1): 0.716065
#> phi (cell 1): 0.998505
#> s (cell 1): 0.607834
#> nu (cell 1): 0.418772
#> theta (batch 1): 0.247981
#> betas:    0.8739  -0.0338   1.3747   0.7091   0.4193   0.3670   0.3151   0.0158  -0.0226  -0.3398  -0.7050  -0.3381
#> 
#> sigma: 0.137828
#> lambda (gene 1): 2.92301
#> epsilon (gene 1): -0.00902563
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 642.284
#> delta (gene 1): 0.860731
#> phi (cell 1): 0.929004
#> s (cell 1): 0.45127
#> nu (cell 1): 0.418597
#> theta (batch 1): 0.229215
#> betas:    0.9735  -0.1197   0.5877   0.8626   0.2660   1.1480   0.0137   0.0526   0.4660  -0.6857   0.2909  -0.3990
#> 
#> sigma: 0.181411
#> lambda (gene 1): 0.49204
#> epsilon (gene 1): -0.0672267
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 736.239
#> delta (gene 1): 0.81494
#> phi (cell 1): 0.838147
#> s (cell 1): 0.691917
#> nu (cell 1): 0.433655
#> theta (batch 1): 0.124765
#> betas:    0.6980  -0.1824   0.3533   0.7783   1.1267   0.6710   0.4285   0.3036   0.2997  -0.0816  -0.0549   0.1575
#> 
#> sigma: 0.119184
#> lambda (gene 1): 0.587967
#> epsilon (gene 1): 0.113062
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 789.12
#> delta (gene 1): 0.62638
#> phi (cell 1): 0.754152
#> s (cell 1): 0.5039
#> nu (cell 1): 0.41775
#> theta (batch 1): 0.290317
#> betas:    1.1367  -0.1017   0.4104   0.7763   0.9829   0.1146   0.3100   0.3159  -0.4016  -0.0668  -0.5451  -0.0993
#> 
#> sigma: 0.151883
#> lambda (gene 1): 0.996564
#> epsilon (gene 1): -0.2336
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 764.992
#> delta (gene 1): 0.556255
#> phi (cell 1): 0.791472
#> s (cell 1): 0.759801
#> nu (cell 1): 0.427894
#> theta (batch 1): 0.207844
#> betas:    0.9301  -0.0453   0.8144   0.1693   0.9476   0.5949  -0.0415   0.2571  -0.4045  -0.2243  -0.1432  -0.9802
#> 
#> sigma: 0.19602
#> lambda (gene 1): 0.653471
#> epsilon (gene 1): -0.38921
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 752.907
#> delta (gene 1): 0.663832
#> phi (cell 1): 0.768661
#> s (cell 1): 0.343058
#> nu (cell 1): 0.414467
#> theta (batch 1): 0.0333121
#> betas:    1.0542  -0.1456   0.5183   0.6710   0.4535   0.6395   0.3105   0.2291  -0.1783  -0.3272   0.2081  -0.2851
#> 
#> sigma: 0.131977
#> lambda (gene 1): 2.76311
#> epsilon (gene 1): -0.190256
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 706.225
#> delta (gene 1): 0.89439
#> phi (cell 1): 0.802838
#> s (cell 1): 0.200895
#> nu (cell 1): 0.41865
#> theta (batch 1): 0.0938644
#> betas:    0.6987  -0.0801   0.1962   0.7464   0.6742   0.4911   0.3091   0.0981   0.0939  -0.1077  -0.4503  -0.2261
#> 
#> sigma: 0.13495
#> lambda (gene 1): 0.32222
#> epsilon (gene 1): 0.145533
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 579.531
#> delta (gene 1): 0.836356
#> phi (cell 1): 0.753836
#> s (cell 1): 0.719317
#> nu (cell 1): 0.420049
#> theta (batch 1): 0.178613
#> betas:    0.9321  -0.0045   0.5617   0.6962   0.0781   0.6847   0.6856  -0.8703   0.2789  -0.4919  -0.8565  -0.4047
#> 
#> sigma: 0.17158
#> lambda (gene 1): 2.34086
#> epsilon (gene 1): 0.157654
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 865.952
#> delta (gene 1): 0.709137
#> phi (cell 1): 0.768282
#> s (cell 1): 0.429128
#> nu (cell 1): 0.411488
#> theta (batch 1): 0.128072
#> betas:    1.3669  -0.0766   0.2578   0.9525   0.0401   0.6303   0.0476  -0.0105  -0.4862   0.0397  -0.3668  -1.0498
#> 
#> sigma: 0.178388
#> lambda (gene 1): 1.73923
#> epsilon (gene 1): -0.339869
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 665.321
#> delta (gene 1): 0.533247
#> phi (cell 1): 0.812006
#> s (cell 1): 0.74675
#> nu (cell 1): 0.407633
#> theta (batch 1): 0.275603
#> betas:    0.1474  -0.0475   1.2450   0.7030   0.6564   0.7774   1.1393  -0.3127   0.4554  -0.3258   0.1509  -0.3420
#> 
#> sigma: 0.205115
#> lambda (gene 1): 0.435315
#> epsilon (gene 1): -0.421719
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 684.219
#> delta (gene 1): 0.585491
#> phi (cell 1): 0.854453
#> s (cell 1): 0.410123
#> nu (cell 1): 0.41918
#> theta (batch 1): 0.0723418
#> betas:    1.3889  -0.2061   0.2612   0.7522   0.5902   0.6279   0.0325   0.0224   0.6579  -0.5205  -0.3709   0.4246
#> 
#> sigma: 0.219141
#> lambda (gene 1): 1.06537
#> epsilon (gene 1): -0.310472
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 620.312
#> delta (gene 1): 0.652942
#> phi (cell 1): 0.788455
#> s (cell 1): 0.383421
#> nu (cell 1): 0.405005
#> theta (batch 1): 0.0442334
#> betas:    0.9674  -0.1091   0.3483   0.7742   0.6912   0.5671   0.1213   0.7174  -0.5248  -0.2373  -0.2000  -0.1022
#> 
#> sigma: 0.121745
#> lambda (gene 1): 1.0621
#> epsilon (gene 1): -0.0661567
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 706.146
#> delta (gene 1): 1.03634
#> phi (cell 1): 0.867682
#> s (cell 1): 0.396822
#> nu (cell 1): 0.426942
#> theta (batch 1): 0.0570178
#> betas:    0.1058   0.1152   1.1949   0.1825   0.4242   1.1427   0.0059   0.4220  -0.2853  -0.3830  -0.7395  -0.3713
#> 
#> sigma: 0.156232
#> lambda (gene 1): 1.7592
#> epsilon (gene 1): 0.30961
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 689.326
#> delta (gene 1): 0.654405
#> phi (cell 1): 0.872174
#> s (cell 1): 1.00664
#> nu (cell 1): 0.408868
#> theta (batch 1): 0.539015
#> betas:    0.8894  -0.1594   0.5708   0.6774   0.1712   0.9938   0.3275   0.2725  -0.0115  -0.1446  -0.0665   0.0827
#> 
#> sigma: 0.137414
#> lambda (gene 1): 0.800429
#> epsilon (gene 1): -0.146409
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 514.18
#> delta (gene 1): 0.897159
#> phi (cell 1): 0.929925
#> s (cell 1): 0.853027
#> nu (cell 1): 0.432527
#> theta (batch 1): 0.0879232
#> betas:    0.7332  -0.0285   0.4955   0.2099   1.2107   0.4112   0.1259  -0.0891   0.3316  -0.3993  -0.2311  -0.6587
#> 
#> sigma: 0.11443
#> lambda (gene 1): 1.25856
#> epsilon (gene 1): -0.0882929
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 748.928
#> delta (gene 1): 0.919764
#> phi (cell 1): 0.977573
#> s (cell 1): 0.428978
#> nu (cell 1): 0.416227
#> theta (batch 1): 0.174272
#> betas:    0.5145   0.0672   0.6124   1.0544   0.0635   1.2600  -0.6511   0.4572  -0.1038  -0.1813  -0.9201  -0.3878
#> 
#> sigma: 0.126439
#> lambda (gene 1): 0.542216
#> epsilon (gene 1): -0.0603845
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 679.869
#> delta (gene 1): 0.646128
#> phi (cell 1): 0.874747
#> s (cell 1): 0.425035
#> nu (cell 1): 0.415302
#> theta (batch 1): 0.121989
#> betas:    1.6201  -0.1608   0.6501   0.5910   0.6989   0.0939   0.1543  -0.0094  -0.1309  -0.1272  -0.3854  -0.6329
#> 
#> sigma: 0.201676
#> lambda (gene 1): 0.2879
#> epsilon (gene 1): -0.300636
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 621.046
#> delta (gene 1): 0.497867
#> phi (cell 1): 0.751481
#> s (cell 1): 0.513358
#> nu (cell 1): 0.420087
#> theta (batch 1): 0.0384665
#> betas:    1.0330  -0.1784   0.5933   0.3899   1.0228   0.7894   0.1347   0.4905   0.0358  -0.1437  -0.0094  -0.0644
#> 
#> sigma: 0.134309
#> lambda (gene 1): 0.52838
#> epsilon (gene 1): -0.508193
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 690.753
#> delta (gene 1): 0.748643
#> phi (cell 1): 0.738383
#> s (cell 1): 0.753093
#> nu (cell 1): 0.4206
#> theta (batch 1): 0.528986
#> betas:    0.8580  -0.1455   1.0174   0.5688   0.8816   0.4146   0.6033   0.0797   0.1308  -0.0473  -0.2377  -0.1060
#> 
#> sigma: 0.158158
#> lambda (gene 1): 1.13848
#> epsilon (gene 1): -0.043559
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 733.594
#> delta (gene 1): 0.730195
#> phi (cell 1): 0.874484
#> s (cell 1): 0.894857
#> nu (cell 1): 0.429859
#> theta (batch 1): 0.441969
#> betas:    1.3152  -0.1143   0.8386   0.5547   0.3016   0.5591   0.4704  -0.3897   0.0857  -0.3231   0.0241  -0.9888
#> 
#> sigma: 0.152666
#> lambda (gene 1): 1.40753
#> epsilon (gene 1): -0.268844
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 792.116
#> delta (gene 1): 0.512198
#> phi (cell 1): 0.785216
#> s (cell 1): 0.390745
#> nu (cell 1): 0.418741
#> theta (batch 1): 0.0944472
#> betas:    0.8906  -0.0923   0.5153   1.0707   0.4739   0.7242   0.1816   0.0686   0.4633  -0.8949   0.2894  -0.4344
#> 
#> sigma: 0.156113
#> lambda (gene 1): 1.1846
#> epsilon (gene 1): -0.472224
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 658.84
#> delta (gene 1): 0.6869
#> phi (cell 1): 0.704342
#> s (cell 1): 0.377142
#> nu (cell 1): 0.422821
#> theta (batch 1): 0.298475
#> betas:    1.5572  -0.2528  -0.0839   1.0177   0.3607   0.9825   0.2350  -0.5461   0.6885  -0.1988  -0.4659   0.5025
#> 
#> sigma: 0.234748
#> lambda (gene 1): 1.16469
#> epsilon (gene 1): -0.186179
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 660.444
#> delta (gene 1): 0.722191
#> phi (cell 1): 0.608161
#> s (cell 1): 0.530289
#> nu (cell 1): 0.422978
#> theta (batch 1): 0.0989868
#> betas:    1.7377  -0.1428   0.3511   0.4263   0.2522   0.7072  -0.0393  -0.1656  -0.1229  -0.5408  -0.4992  -0.3069
#> 
#> sigma: 0.189412
#> lambda (gene 1): 2.42826
#> epsilon (gene 1): 0.00932264
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 784.899
#> delta (gene 1): 0.638058
#> phi (cell 1): 0.672792
#> s (cell 1): 1.05984
#> nu (cell 1): 0.428932
#> theta (batch 1): 0.282446
#> betas:    1.4509  -0.1464   0.0848   0.2829   0.9541   0.4693  -0.1740   0.2040  -0.0633  -0.1493  -0.6368   0.0741
#> 
#> sigma: 0.243199
#> lambda (gene 1): 0.534284
#> epsilon (gene 1): -0.310912
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 849.975
#> delta (gene 1): 0.974504
#> phi (cell 1): 0.82948
#> s (cell 1): 1.81158
#> nu (cell 1): 0.418735
#> theta (batch 1): 0.552659
#> betas:    0.2850  -0.0502   0.8552   0.5037   0.7930   0.9057   0.3884   0.0923   0.2433  -0.1149  -0.3826   0.0474
#> 
#> sigma: 0.121651
#> lambda (gene 1): 2.15122
#> epsilon (gene 1): 0.222422
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 669.981
#> delta (gene 1): 0.899051
#> phi (cell 1): 0.644499
#> s (cell 1): 0.324596
#> nu (cell 1): 0.433049
#> theta (batch 1): 0.21046
#> betas:    0.3646  -0.0065   0.4870   1.3257   0.4658   0.8776   0.3984   0.0990   0.1101  -0.4609   0.1099  -0.6271
#> 
#> sigma: 0.163458
#> lambda (gene 1): 1.72561
#> epsilon (gene 1): 0.021923
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 632.799
#> delta (gene 1): 0.710377
#> phi (cell 1): 0.687342
#> s (cell 1): 0.788505
#> nu (cell 1): 0.406794
#> theta (batch 1): 0.394758
#> betas:    1.0426  -0.1216   0.6318   0.4362   0.6412   0.5706   0.2307   0.2993   0.1360  -0.2063  -0.4068  -0.3031
#> 
#> sigma: 0.175937
#> lambda (gene 1): 1.85925
#> epsilon (gene 1): -0.132636
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 700.004
#> delta (gene 1): 0.591879
#> phi (cell 1): 0.748979
#> s (cell 1): 0.843379
#> nu (cell 1): 0.414154
#> theta (batch 1): 0.200639
#> betas:    0.7092  -0.0477   0.6420   1.1952   0.3468   0.5085   0.7237  -0.1897  -0.0758  -0.5295   0.1355  -0.6643
#> 
#> sigma: 0.151901
#> lambda (gene 1): 1.18713
#> epsilon (gene 1): -0.23612
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 813.728
#> delta (gene 1): 0.679481
#> phi (cell 1): 0.637644
#> s (cell 1): 0.3697
#> nu (cell 1): 0.426281
#> theta (batch 1): 0.303455
#> betas:    0.5842  -0.0758   0.5546   0.7857   1.1411   0.3559   0.2887   0.4099   0.0033  -0.3290  -0.0772  -0.0153
#> 
#> sigma: 0.174016
#> lambda (gene 1): 0.415645
#> epsilon (gene 1): -0.175171
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 649.001
#> delta (gene 1): 0.553093
#> phi (cell 1): 0.729102
#> s (cell 1): 0.302233
#> nu (cell 1): 0.414201
#> theta (batch 1): 0.360524
#> betas:    0.6016   0.0110   0.7022   1.3925  -0.1489   0.5013   0.3967   0.1610  -0.1988  -0.3731  -0.3164  -0.6348
#> 
#> sigma: 0.157847
#> lambda (gene 1): 0.449813
#> epsilon (gene 1): -0.389054
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 869.834
#> delta (gene 1): 0.630403
#> phi (cell 1): 0.739085
#> s (cell 1): 2.21321
#> nu (cell 1): 0.413399
#> theta (batch 1): 0.200454
#> betas:    1.1583  -0.0906   0.7445   0.1549   1.1991   0.1957   0.0029   0.1771   0.0276  -0.4021  -0.5361  -0.0807
#> 
#> sigma: 0.140339
#> lambda (gene 1): 0.980166
#> epsilon (gene 1): -0.240571
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 742.184
#> delta (gene 1): 0.974496
#> phi (cell 1): 0.686067
#> s (cell 1): 0.422934
#> nu (cell 1): 0.42128
#> theta (batch 1): 0.430449
#> betas:    1.3439  -0.2595   0.9960   0.4628   0.6197   0.0865   0.7356   0.3595  -0.2661   0.3290  -0.2005   0.1640
#> 
#> sigma: 0.127116
#> lambda (gene 1): 0.530411
#> epsilon (gene 1): 0.205575
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 759.342
#> delta (gene 1): 0.635912
#> phi (cell 1): 0.778409
#> s (cell 1): 0.763026
#> nu (cell 1): 0.420933
#> theta (batch 1): 0.258067
#> betas:    0.6261   0.0397   0.5104   0.5226   0.9615   0.4259  -0.1477   0.1476  -0.0981  -0.4432  -0.7174  -0.4477
#> 
#> sigma: 0.128151
#> lambda (gene 1): 0.567581
#> epsilon (gene 1): -0.187938
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 606.193
#> delta (gene 1): 0.706889
#> phi (cell 1): 0.823864
#> s (cell 1): 0.383518
#> nu (cell 1): 0.425728
#> theta (batch 1): 0.0197522
#> betas:    1.0062  -0.0735   0.6558   0.1894   1.1935   0.0510   0.4494  -0.0235  -0.0794  -0.4257  -0.0925  -0.4887
#> 
#> sigma: 0.161071
#> lambda (gene 1): 0.684901
#> epsilon (gene 1): -0.186037
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 741.921
#> delta (gene 1): 0.704227
#> phi (cell 1): 0.846704
#> s (cell 1): 0.454802
#> nu (cell 1): 0.423782
#> theta (batch 1): 0.0765913
#> betas:    1.4007  -0.1783   0.8527  -0.0116   0.3973   1.1093  -0.2552   0.0802  -0.0658  -0.0762  -0.4731   0.3327
#> 
#> sigma: 0.163799
#> lambda (gene 1): 0.658274
#> epsilon (gene 1): -0.194018
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 541.179
#> delta (gene 1): 0.581579
#> phi (cell 1): 0.736727
#> s (cell 1): 0.729197
#> nu (cell 1): 0.425254
#> theta (batch 1): 0.165758
#> betas:    1.5463  -0.0821   0.0311   0.1990   0.8733   0.1307   0.2544  -0.2164  -0.0953  -0.5267  -0.7939  -0.1434
#> 
#> sigma: 0.155221
#> lambda (gene 1): 0.544583
#> epsilon (gene 1): -0.272491
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 631.167
#> delta (gene 1): 0.689301
#> phi (cell 1): 0.685111
#> s (cell 1): 0.870543
#> nu (cell 1): 0.445533
#> theta (batch 1): 0.104218
#> betas:    0.6243   0.0163   0.6634   0.7741   0.9684  -0.0499   0.5415   0.0741  -0.0441  -0.7716  -0.1071  -0.3007
#> 
#> sigma: 0.176149
#> lambda (gene 1): 1.1188
#> epsilon (gene 1): -0.173101
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 726.41
#> delta (gene 1): 0.605395
#> phi (cell 1): 0.85095
#> s (cell 1): 0.274382
#> nu (cell 1): 0.42355
#> theta (batch 1): 0.206744
#> betas:    1.0292  -0.1008   0.7420   0.3847   0.6148   0.9736   0.2227  -0.3440   0.1586  -0.0004  -0.6975  -0.0117
#> 
#> sigma: 0.126348
#> lambda (gene 1): 1.19659
#> epsilon (gene 1): -0.365038
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 575.195
#> delta (gene 1): 0.798825
#> phi (cell 1): 0.892377
#> s (cell 1): 0.607834
#> nu (cell 1): 0.41318
#> theta (batch 1): 0.0741937
#> betas:    0.5199   0.0249   0.1944   0.7017   0.7919   0.3381   0.8215  -0.3302  -0.4077  -0.0911  -0.4762  -0.2683
#> 
#> sigma: 0.129731
#> lambda (gene 1): 1.29381
#> epsilon (gene 1): -0.0732441
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 673.767
#> delta (gene 1): 0.588963
#> phi (cell 1): 0.968266
#> s (cell 1): 0.237081
#> nu (cell 1): 0.429983
#> theta (batch 1): 0.0922194
#> betas:    0.9574  -0.1689   0.4688   0.6338   1.1743   0.3403   0.5811   0.0962   0.2040   0.1220  -0.4781   0.3220
#> 
#> sigma: 0.197855
#> lambda (gene 1): 0.940407
#> epsilon (gene 1): -0.387146
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 722.402
#> delta (gene 1): 0.957825
#> phi (cell 1): 1.01246
#> s (cell 1): 0.770091
#> nu (cell 1): 0.429781
#> theta (batch 1): 0.112167
#> betas:    1.2168  -0.1091   0.3803   0.6758   0.1718   1.0101  -0.2222   0.3850  -0.0674  -0.3930  -0.3396  -0.4608
#> 
#> sigma: 0.14364
#> lambda (gene 1): 1.46914
#> epsilon (gene 1): 0.218061
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 708.432
#> delta (gene 1): 0.590575
#> phi (cell 1): 0.923286
#> s (cell 1): 0.337799
#> nu (cell 1): 0.412218
#> theta (batch 1): 0.246321
#> betas:    0.7677  -0.0647   0.4948   0.6884   0.5075   0.4927   0.3376   0.1273  -0.0354  -0.3279  -0.1712  -0.3247
#> 
#> sigma: 0.147999
#> lambda (gene 1): 0.276888
#> epsilon (gene 1): -0.329285
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 648.453
#> delta (gene 1): 0.546009
#> phi (cell 1): 0.915697
#> s (cell 1): 0.362136
#> nu (cell 1): 0.406417
#> theta (batch 1): 0.186708
#> betas:    1.1916  -0.1180   0.8815   0.4281   0.7161   0.3798   0.3316  -0.0649   0.3845  -0.4175  -0.5531   0.1230
#> 
#> sigma: 0.208416
#> lambda (gene 1): 0.45216
#> epsilon (gene 1): -0.458091
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 623.939
#> delta (gene 1): 0.491528
#> phi (cell 1): 0.801659
#> s (cell 1): 0.36444
#> nu (cell 1): 0.434566
#> theta (batch 1): 0.506554
#> betas:    0.9925  -0.1941   0.3579   1.1160   0.5767   0.1507   0.9663   0.2682  -0.1868   0.1611  -0.3213   0.1946
#> 
#> sigma: 0.18747
#> lambda (gene 1): 1.78363
#> epsilon (gene 1): -0.352106
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 860.943
#> delta (gene 1): 0.558085
#> phi (cell 1): 0.748891
#> s (cell 1): 0.61519
#> nu (cell 1): 0.417021
#> theta (batch 1): 0.111974
#> betas:    0.7072  -0.0563   0.3322   0.9386   1.1933   0.5057   0.3335  -0.0351  -0.0100  -0.2757  -0.6137   0.0472
#> 
#> sigma: 0.185795
#> lambda (gene 1): 1.03967
#> epsilon (gene 1): -0.186405
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 808.397
#> delta (gene 1): 0.650197
#> phi (cell 1): 0.837732
#> s (cell 1): 0.51956
#> nu (cell 1): 0.419576
#> theta (batch 1): 0.118997
#> betas:    0.7049  -0.0961   0.4425   0.8247   0.6791   0.5550   0.2262   0.7251  -0.4109   0.1973  -0.4212  -0.0559
#> 
#> sigma: 0.124955
#> lambda (gene 1): 1.018
#> epsilon (gene 1): -0.250429
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 5980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 761.614
#> delta (gene 1): 0.741122
#> phi (cell 1): 0.944922
#> s (cell 1): 0.658812
#> nu (cell 1): 0.408842
#> theta (batch 1): 0.243444
#> betas:    1.4114  -0.2119   0.4995   0.7028   0.4001   0.3992   0.4209  -0.0599   0.6278  -0.6686   0.0195  -0.2329
#> 
#> sigma: 0.144037
#> lambda (gene 1): 0.799304
#> epsilon (gene 1): 0.00728885
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 567.201
#> delta (gene 1): 0.645593
#> phi (cell 1): 0.900382
#> s (cell 1): 0.641074
#> nu (cell 1): 0.422267
#> theta (batch 1): 0.379915
#> betas:    1.0112  -0.1303   0.7913   0.4922   0.7443   0.7170   0.2917  -0.3599   0.6346  -0.2734  -0.1856  -0.4108
#> 
#> sigma: 0.123545
#> lambda (gene 1): 0.550224
#> epsilon (gene 1): -0.396111
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 667.429
#> delta (gene 1): 0.795956
#> phi (cell 1): 0.76537
#> s (cell 1): 0.315966
#> nu (cell 1): 0.411871
#> theta (batch 1): 0.109966
#> betas:    1.1345  -0.0335   0.4744   0.5480   0.5750   0.4938   0.2933   0.0075  -0.6036  -0.0386  -0.4650  -0.7799
#> 
#> sigma: 0.158572
#> lambda (gene 1): 0.342369
#> epsilon (gene 1): -0.140673
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 821.58
#> delta (gene 1): 0.857167
#> phi (cell 1): 0.760005
#> s (cell 1): 0.650712
#> nu (cell 1): 0.427668
#> theta (batch 1): 0.605447
#> betas:    0.5045  -0.0242   0.6345   0.8633   0.2552   1.0877   0.6508  -0.4118   0.0057   0.2573  -0.8929  -0.3171
#> 
#> sigma: 0.224959
#> lambda (gene 1): 0.798139
#> epsilon (gene 1): 0.077927
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1028
#> delta (gene 1): 0.53445
#> phi (cell 1): 0.757388
#> s (cell 1): 0.460778
#> nu (cell 1): 0.42034
#> theta (batch 1): 0.201378
#> betas:    0.9689  -0.1879   0.6146   0.7864   0.3130   1.1009   0.1385   0.0851   0.8125  -0.8801   0.3510  -0.0235
#> 
#> sigma: 0.160047
#> lambda (gene 1): 1.02319
#> epsilon (gene 1): -0.282926
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 577.498
#> delta (gene 1): 0.78611
#> phi (cell 1): 0.889186
#> s (cell 1): 0.187683
#> nu (cell 1): 0.403067
#> theta (batch 1): 0.255941
#> betas:    0.3538  -0.0511   0.1872   1.4069   0.4264   0.8790   0.5574   0.0200  -0.0340   0.1866  -0.5734   0.1131
#> 
#> sigma: 0.142408
#> lambda (gene 1): 0.737676
#> epsilon (gene 1): -0.0395912
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 808.134
#> delta (gene 1): 0.636521
#> phi (cell 1): 0.888421
#> s (cell 1): 0.908846
#> nu (cell 1): 0.413146
#> theta (batch 1): 0.141961
#> betas:    1.0544  -0.1218   0.5254   0.3034   0.4721   0.6787   0.5849  -0.4741   0.3995  -0.1570  -0.5363   0.0866
#> 
#> sigma: 0.161249
#> lambda (gene 1): 0.669133
#> epsilon (gene 1): -0.336755
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 612.865
#> delta (gene 1): 0.60523
#> phi (cell 1): 0.905674
#> s (cell 1): 0.670223
#> nu (cell 1): 0.416218
#> theta (batch 1): 0.346673
#> betas:    1.5211  -0.1243  -0.0092   0.3905   0.7411   0.2481   0.1488   0.2468  -0.6566   0.1286  -0.7305  -0.2310
#> 
#> sigma: 0.141214
#> lambda (gene 1): 0.733831
#> epsilon (gene 1): -0.3315
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 805.33
#> delta (gene 1): 0.586632
#> phi (cell 1): 0.853159
#> s (cell 1): 0.407338
#> nu (cell 1): 0.419637
#> theta (batch 1): 0.318832
#> betas:    0.5308  -0.1105   0.7153   1.0476   0.6124   0.3589   0.7953   0.0781   0.1741   0.0235   0.1198  -0.4186
#> 
#> sigma: 0.161681
#> lambda (gene 1): 0.632806
#> epsilon (gene 1): -0.479342
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 745.58
#> delta (gene 1): 0.864718
#> phi (cell 1): 0.788397
#> s (cell 1): 0.413644
#> nu (cell 1): 0.421538
#> theta (batch 1): 0.230418
#> betas:    1.2933  -0.0864   0.6082   0.3653   0.6646   0.5596  -0.0717   0.0730  -0.3365   0.2484  -0.7979  -0.3809
#> 
#> sigma: 0.172753
#> lambda (gene 1): 1.22117
#> epsilon (gene 1): -0.187939
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 605.974
#> delta (gene 1): 0.690136
#> phi (cell 1): 0.803171
#> s (cell 1): 0.468557
#> nu (cell 1): 0.399095
#> theta (batch 1): 0.262524
#> betas:    0.1524  -0.0273   0.5171   1.4018   0.6748   0.6922   0.2471   0.7836  -0.3257  -0.1693  -0.0207  -0.2317
#> 
#> sigma: 0.159376
#> lambda (gene 1): 0.944834
#> epsilon (gene 1): -0.0546723
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 741.817
#> delta (gene 1): 0.756173
#> phi (cell 1): 0.726567
#> s (cell 1): 0.203266
#> nu (cell 1): 0.42041
#> theta (batch 1): 0.49315
#> betas:    0.9052   0.0055   0.6650   0.5354   0.3699   1.1168  -0.2424   0.0800  -0.2129  -0.7600  -0.1694  -0.7343
#> 
#> sigma: 0.150246
#> lambda (gene 1): 3.52862
#> epsilon (gene 1): 0.00508362
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 633.364
#> delta (gene 1): 0.574921
#> phi (cell 1): 0.801258
#> s (cell 1): 0.791813
#> nu (cell 1): 0.41788
#> theta (batch 1): 0.215465
#> betas:    1.5623  -0.1512   0.4589   0.4733   0.4328   0.6251   0.1233   0.2224  -0.4130   0.0935  -1.1994   0.5275
#> 
#> sigma: 0.169314
#> lambda (gene 1): 0.712861
#> epsilon (gene 1): -0.227173
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 625.405
#> delta (gene 1): 0.746857
#> phi (cell 1): 0.715957
#> s (cell 1): 0.999833
#> nu (cell 1): 0.425746
#> theta (batch 1): 0.582819
#> betas:    0.6743  -0.0561   0.6781   0.7979   0.9728   0.1130   0.3283   0.6717  -0.2691  -0.2851  -0.4700   0.0720
#> 
#> sigma: 0.190212
#> lambda (gene 1): 1.19117
#> epsilon (gene 1): 0.0495308
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 637.37
#> delta (gene 1): 0.498261
#> phi (cell 1): 0.714509
#> s (cell 1): 1.09921
#> nu (cell 1): 0.406659
#> theta (batch 1): 0.500018
#> betas:    0.7612  -0.0636   1.1496   0.8858   0.5419   0.6309   0.5484  -0.0947  -0.1452  -0.3593  -0.0924  -0.2797
#> 
#> sigma: 0.112572
#> lambda (gene 1): 0.665103
#> epsilon (gene 1): -0.431336
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 655.638
#> delta (gene 1): 0.636496
#> phi (cell 1): 0.705434
#> s (cell 1): 0.535604
#> nu (cell 1): 0.429755
#> theta (batch 1): 0.449112
#> betas:    0.8854  -0.1391   0.6939   0.7503   0.4828   0.8954   0.0940   0.2556   0.0893  -0.0920  -0.4577   0.2598
#> 
#> sigma: 0.145818
#> lambda (gene 1): 0.720821
#> epsilon (gene 1): -0.164732
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 763.412
#> delta (gene 1): 0.601335
#> phi (cell 1): 0.67305
#> s (cell 1): 0.556822
#> nu (cell 1): 0.421198
#> theta (batch 1): 0.380346
#> betas:    1.3364  -0.1701   0.7805   0.1490   0.6443   0.5358   0.0943   0.1378  -0.0742   0.0239  -0.6023   0.0032
#> 
#> sigma: 0.162447
#> lambda (gene 1): 1.57667
#> epsilon (gene 1): -0.264226
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 696.932
#> delta (gene 1): 0.596211
#> phi (cell 1): 0.621069
#> s (cell 1): 0.393419
#> nu (cell 1): 0.420013
#> theta (batch 1): 0.0214306
#> betas:    0.5110   0.0379  -0.0944   1.3924   0.7657   0.4637   0.1544  -0.0775   0.3145  -0.9135  -0.1626  -0.7844
#> 
#> sigma: 0.144646
#> lambda (gene 1): 3.90977
#> epsilon (gene 1): -0.212767
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 639.62
#> delta (gene 1): 0.867138
#> phi (cell 1): 0.61868
#> s (cell 1): 0.35213
#> nu (cell 1): 0.403265
#> theta (batch 1): 0.0215306
#> betas:    1.3357  -0.2143   0.9419   0.5348   0.1553   0.6946   0.5532   0.0168   0.0822  -0.2362   0.2130  -0.2467
#> 
#> sigma: 0.163372
#> lambda (gene 1): 0.199649
#> epsilon (gene 1): 0.0147993
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 646.08
#> delta (gene 1): 0.856102
#> phi (cell 1): 0.631939
#> s (cell 1): 0.445886
#> nu (cell 1): 0.43528
#> theta (batch 1): 0.0164159
#> betas:    0.8127  -0.0766   0.7636   0.8663   0.3972   0.4742   0.3581   0.0267   0.0709  -0.0976  -0.8343   0.2927
#> 
#> sigma: 0.137262
#> lambda (gene 1): 1.6013
#> epsilon (gene 1): 0.124513
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 642.096
#> delta (gene 1): 0.653373
#> phi (cell 1): 0.692019
#> s (cell 1): 1.26894
#> nu (cell 1): 0.409592
#> theta (batch 1): 0.246091
#> betas:    0.9717  -0.2424   0.7545   1.0537   0.3601   1.1134   0.1391   0.3497   0.1151   0.3896  -0.0431   0.2333
#> 
#> sigma: 0.154455
#> lambda (gene 1): 0.620216
#> epsilon (gene 1): -0.363956
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 802.835
#> delta (gene 1): 0.598025
#> phi (cell 1): 0.572785
#> s (cell 1): 0.419846
#> nu (cell 1): 0.419594
#> theta (batch 1): 0.295364
#> betas:    0.8245  -0.0822  -0.0257   1.1138   0.6533   0.6764  -0.0097   0.2979  -0.0312   0.1026  -0.6618  -0.3031
#> 
#> sigma: 0.116835
#> lambda (gene 1): 2.09974
#> epsilon (gene 1): -0.345031
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 801.393
#> delta (gene 1): 0.842386
#> phi (cell 1): 0.643875
#> s (cell 1): 0.589588
#> nu (cell 1): 0.400615
#> theta (batch 1): 0.289369
#> betas:    1.3068  -0.2310   0.3098   0.1629   0.7935   0.7841   0.1668   0.3120   0.2679  -0.1100  -0.6375   0.5737
#> 
#> sigma: 0.163775
#> lambda (gene 1): 1.2281
#> epsilon (gene 1): 0.239194
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 737.406
#> delta (gene 1): 0.420033
#> phi (cell 1): 0.757446
#> s (cell 1): 0.381574
#> nu (cell 1): 0.416425
#> theta (batch 1): 0.201068
#> betas:    0.6292   0.0196   0.8433   0.6398   1.0056   0.2124   0.3026   0.1064  -0.2956  -0.1020  -0.6219  -0.6962
#> 
#> sigma: 0.142935
#> lambda (gene 1): 0.476828
#> epsilon (gene 1): -0.690153
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 560.258
#> delta (gene 1): 0.542979
#> phi (cell 1): 0.831912
#> s (cell 1): 0.362638
#> nu (cell 1): 0.411801
#> theta (batch 1): 0.183654
#> betas:    1.0028  -0.0666   0.6578   0.7499   0.8142   0.0198   0.4539   0.1287  -0.0906  -0.2349  -0.4311  -0.5728
#> 
#> sigma: 0.133943
#> lambda (gene 1): 0.766955
#> epsilon (gene 1): -0.424994
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 641.093
#> delta (gene 1): 0.702826
#> phi (cell 1): 0.799171
#> s (cell 1): 0.328358
#> nu (cell 1): 0.411801
#> theta (batch 1): 0.404111
#> betas:    0.9487  -0.1799   0.9349   0.7922   0.3365   0.5184   0.9720  -0.0169  -0.0779   0.1154  -0.1085   0.0061
#> 
#> sigma: 0.215413
#> lambda (gene 1): 0.727564
#> epsilon (gene 1): -0.119002
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 474.651
#> delta (gene 1): 0.738382
#> phi (cell 1): 0.735378
#> s (cell 1): 0.273365
#> nu (cell 1): 0.422465
#> theta (batch 1): 0.148014
#> betas:    1.5604  -0.0637   0.8639   0.4131   0.2295   0.6278  -0.3487  -0.1348   0.0347  -0.7453  -0.8730  -0.3185
#> 
#> sigma: 0.184985
#> lambda (gene 1): 0.458893
#> epsilon (gene 1): 0.0768837
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 564.985
#> delta (gene 1): 0.515131
#> phi (cell 1): 0.787898
#> s (cell 1): 0.476322
#> nu (cell 1): 0.421198
#> theta (batch 1): 0.248274
#> betas:    0.7516  -0.0990   0.4924   0.5386   0.9600   0.0852   0.9422   0.0215   0.1611  -0.4505  -0.3262   0.1508
#> 
#> sigma: 0.130526
#> lambda (gene 1): 0.621143
#> epsilon (gene 1): -0.244828
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 541.511
#> delta (gene 1): 0.574495
#> phi (cell 1): 0.75111
#> s (cell 1): 0.291303
#> nu (cell 1): 0.431984
#> theta (batch 1): 0.308053
#> betas:    0.3384  -0.0172   1.3058   0.2709   1.0536   0.8019   0.2993   0.0412   0.3405  -0.4771  -0.1939  -0.3706
#> 
#> sigma: 0.145893
#> lambda (gene 1): 0.565279
#> epsilon (gene 1): -0.289772
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 567.21
#> delta (gene 1): 0.60833
#> phi (cell 1): 0.699232
#> s (cell 1): 0.621095
#> nu (cell 1): 0.421568
#> theta (batch 1): 0.203092
#> betas:    1.2487  -0.1999   0.2833   0.6512   0.8018   0.3626  -0.0277   0.3293   0.2039  -0.2160   0.0431   0.1657
#> 
#> sigma: 0.218606
#> lambda (gene 1): 0.925012
#> epsilon (gene 1): -0.483572
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 556.629
#> delta (gene 1): 0.891652
#> phi (cell 1): 0.810876
#> s (cell 1): 0.917242
#> nu (cell 1): 0.421568
#> theta (batch 1): 0.602066
#> betas:    0.7471  -0.0796   0.4155   0.8190   1.2156   0.0889   0.7722  -0.2561   0.3487  -0.4663  -0.0529  -0.3078
#> 
#> sigma: 0.162822
#> lambda (gene 1): 0.612051
#> epsilon (gene 1): 0.0634404
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 713.491
#> delta (gene 1): 0.836813
#> phi (cell 1): 0.801744
#> s (cell 1): 0.374593
#> nu (cell 1): 0.413091
#> theta (batch 1): 0.23293
#> betas:    0.7093   0.0828   0.1034   0.8971   0.7800   0.4193   0.1295  -0.1797  -0.3186  -0.2215  -0.7838  -0.9211
#> 
#> sigma: 0.125763
#> lambda (gene 1): 1.83077
#> epsilon (gene 1): -0.112136
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 731.682
#> delta (gene 1): 0.811891
#> phi (cell 1): 0.771967
#> s (cell 1): 0.389271
#> nu (cell 1): 0.431469
#> theta (batch 1): 0.140993
#> betas:    0.6317  -0.1131   0.5662   0.6873   0.5599   0.8014   0.4217  -0.1858   0.6740  -0.2352  -0.5700   0.5665
#> 
#> sigma: 0.242423
#> lambda (gene 1): 1.26459
#> epsilon (gene 1): 0.00101828
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 633.893
#> delta (gene 1): 0.602273
#> phi (cell 1): 0.756139
#> s (cell 1): 0.411558
#> nu (cell 1): 0.420705
#> theta (batch 1): 0.330641
#> betas:    0.6993  -0.0887   0.7843   0.8320   0.7615   0.4730   0.4662  -0.0322   0.5520  -0.5883  -0.3674   0.2216
#> 
#> sigma: 0.16122
#> lambda (gene 1): 1.23996
#> epsilon (gene 1): -0.183611
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 712.828
#> delta (gene 1): 0.721411
#> phi (cell 1): 0.775922
#> s (cell 1): 0.62587
#> nu (cell 1): 0.423038
#> theta (batch 1): 0.157318
#> betas:    0.1538   0.0009   0.1912   1.1812   0.6920   0.3668   0.6139   0.3559   0.3358  -0.8041  -0.0552  -0.1243
#> 
#> sigma: 0.151739
#> lambda (gene 1): 0.806303
#> epsilon (gene 1): 0.0773466
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 625.699
#> delta (gene 1): 0.642202
#> phi (cell 1): 0.69879
#> s (cell 1): 0.435981
#> nu (cell 1): 0.4229
#> theta (batch 1): 0.21652
#> betas:    0.4109  -0.0946   0.8080   0.9542   0.5462   1.1332   0.1398   0.3069   0.2554  -0.2155  -0.0091  -0.0191
#> 
#> sigma: 0.164138
#> lambda (gene 1): 2.2655
#> epsilon (gene 1): -0.235681
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 744.579
#> delta (gene 1): 0.728457
#> phi (cell 1): 0.765395
#> s (cell 1): 0.583926
#> nu (cell 1): 0.429926
#> theta (batch 1): 0.306163
#> betas:    0.9068  -0.1862   0.5441   0.5399   0.7100   0.3955   0.5113   0.3037   0.2481  -0.1658  -0.1786   0.4498
#> 
#> sigma: 0.187684
#> lambda (gene 1): 1.8879
#> epsilon (gene 1): -0.0582843
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 760.768
#> delta (gene 1): 0.633495
#> phi (cell 1): 0.770937
#> s (cell 1): 1.50842
#> nu (cell 1): 0.416575
#> theta (batch 1): 0.516764
#> betas:    1.2753  -0.1636   1.1082   0.1102   0.6881   0.4210   0.2242  -0.0583   0.0738   0.1067  -0.4894   0.2325
#> 
#> sigma: 0.18577
#> lambda (gene 1): 1.13558
#> epsilon (gene 1): -0.488686
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 836.321
#> delta (gene 1): 0.879561
#> phi (cell 1): 0.834556
#> s (cell 1): 0.639948
#> nu (cell 1): 0.41247
#> theta (batch 1): 0.191959
#> betas:    0.9732  -0.2164   0.1348   0.9081   1.0119   0.5785   0.1635   0.4903   0.0276   0.0054   0.3241  -0.3172
#> 
#> sigma: 0.239409
#> lambda (gene 1): 1.34073
#> epsilon (gene 1): 0.0496316
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 603.441
#> delta (gene 1): 0.671117
#> phi (cell 1): 0.871841
#> s (cell 1): 0.45518
#> nu (cell 1): 0.416517
#> theta (batch 1): 0.497726
#> betas:    0.6994  -0.1523   0.4875   1.1644   0.3736   1.1527   0.1208   0.1530   0.7552  -0.5594   0.0067   0.0249
#> 
#> sigma: 0.187705
#> lambda (gene 1): 2.95868
#> epsilon (gene 1): -0.0699974
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 799.816
#> delta (gene 1): 0.810702
#> phi (cell 1): 0.873138
#> s (cell 1): 0.637573
#> nu (cell 1): 0.424301
#> theta (batch 1): 0.165798
#> betas:    0.8581   0.0147   0.6847   0.5681   0.5893   0.3088   0.4425  -0.2015  -0.5173   0.2819  -1.3815  -0.4037
#> 
#> sigma: 0.178993
#> lambda (gene 1): 1.76744
#> epsilon (gene 1): 0.0871319
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 592.784
#> delta (gene 1): 0.879866
#> phi (cell 1): 0.871744
#> s (cell 1): 0.889182
#> nu (cell 1): 0.40386
#> theta (batch 1): 0.446217
#> betas:    0.8997  -0.1959   0.1637   1.1872   0.7486   0.6039   0.5109   0.5315  -0.2591   0.3721   0.0959  -0.1848
#> 
#> sigma: 0.160677
#> lambda (gene 1): 2.13234
#> epsilon (gene 1): -0.0774873
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 636.007
#> delta (gene 1): 0.649897
#> phi (cell 1): 0.802849
#> s (cell 1): 0.512102
#> nu (cell 1): 0.407199
#> theta (batch 1): 0.27542
#> betas:    0.5018  -0.1462   0.8445   1.1335   0.8462   0.8069  -0.0036   0.6852   0.0515   0.1573   0.2006  -0.5178
#> 
#> sigma: 0.236552
#> lambda (gene 1): 2.32402
#> epsilon (gene 1): -0.298649
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 671.407
#> delta (gene 1): 0.893216
#> phi (cell 1): 0.899896
#> s (cell 1): 0.821172
#> nu (cell 1): 0.409304
#> theta (batch 1): 0.316368
#> betas:    0.8636  -0.0349   0.3235   0.6719   0.4868   0.6389   0.1127   0.0939  -0.0292  -0.4914  -0.0622  -0.9221
#> 
#> sigma: 0.155888
#> lambda (gene 1): 0.90922
#> epsilon (gene 1): 0.0546672
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 604.499
#> delta (gene 1): 0.825471
#> phi (cell 1): 0.821328
#> s (cell 1): 0.239456
#> nu (cell 1): 0.418209
#> theta (batch 1): 0.0641371
#> betas:    0.8993   0.0488   0.3329   1.0710   0.2139   0.3666   0.1630  -0.4185   0.2406  -0.5868  -1.0646  -0.2695
#> 
#> sigma: 0.129058
#> lambda (gene 1): 0.967389
#> epsilon (gene 1): -0.0158922
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 590.212
#> delta (gene 1): 0.533659
#> phi (cell 1): 0.884484
#> s (cell 1): 0.412331
#> nu (cell 1): 0.429483
#> theta (batch 1): 0.219873
#> betas:    0.5541  -0.0357   0.9370   0.7208   0.7421   0.8524   0.1707  -0.0307   0.1243  -0.0071  -0.3390  -0.4877
#> 
#> sigma: 0.164552
#> lambda (gene 1): 0.195522
#> epsilon (gene 1): -0.618954
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 651.087
#> delta (gene 1): 0.612864
#> phi (cell 1): 0.90902
#> s (cell 1): 0.629338
#> nu (cell 1): 0.433216
#> theta (batch 1): 0.149725
#> betas:    0.1673  -0.0903   0.7686   1.5426   0.6581   0.7577   0.4849   0.4678   0.5798  -0.6788   0.5226  -0.0583
#> 
#> sigma: 0.193244
#> lambda (gene 1): 1.22012
#> epsilon (gene 1): -0.25577
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 623.818
#> delta (gene 1): 0.6083
#> phi (cell 1): 0.894677
#> s (cell 1): 0.631236
#> nu (cell 1): 0.424556
#> theta (batch 1): 0.255964
#> betas:   -0.1136   0.0391   0.4420   0.9829   1.3901   0.5816   0.3567   0.3450   0.2233  -0.2040  -0.4777  -0.2148
#> 
#> sigma: 0.177949
#> lambda (gene 1): 1.95058
#> epsilon (gene 1): -0.245315
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 814.194
#> delta (gene 1): 0.950134
#> phi (cell 1): 0.893357
#> s (cell 1): 0.411987
#> nu (cell 1): 0.433113
#> theta (batch 1): 0.429654
#> betas:    1.2666  -0.1693   0.1919   0.2909   0.9253   0.3977   0.5216  -0.2073   0.4298  -0.3491  -0.3209  -0.1846
#> 
#> sigma: 0.173139
#> lambda (gene 1): 0.859245
#> epsilon (gene 1): 0.200423
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 699.954
#> delta (gene 1): 1.07382
#> phi (cell 1): 0.945116
#> s (cell 1): 0.684795
#> nu (cell 1): 0.420329
#> theta (batch 1): 0.343157
#> betas:    0.5455  -0.0572   0.5397   0.7499   0.6648   0.9517   0.4814  -0.2667   0.2530  -0.3527   0.1967  -0.6923
#> 
#> sigma: 0.177507
#> lambda (gene 1): 0.908049
#> epsilon (gene 1): 0.182092
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 6980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 876.881
#> delta (gene 1): 0.693142
#> phi (cell 1): 0.809105
#> s (cell 1): 0.492571
#> nu (cell 1): 0.43225
#> theta (batch 1): 0.0980936
#> betas:    0.5393  -0.1282   0.3249   0.9011   1.1775   0.3245   0.7321   0.2682   0.0917  -0.0888   0.0402   0.0226
#> 
#> sigma: 0.213255
#> lambda (gene 1): 1.27204
#> epsilon (gene 1): -0.143627
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 705.784
#> delta (gene 1): 0.743656
#> phi (cell 1): 0.786968
#> s (cell 1): 0.418177
#> nu (cell 1): 0.421419
#> theta (batch 1): 0.221335
#> betas:    1.0641  -0.2064   1.0213   0.3854   0.7892   0.5625   0.2633   0.5695   0.0843  -0.1808   0.1867  -0.3116
#> 
#> sigma: 0.182815
#> lambda (gene 1): 0.28821
#> epsilon (gene 1): -0.0438864
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 845.025
#> delta (gene 1): 0.921039
#> phi (cell 1): 0.702054
#> s (cell 1): 0.316669
#> nu (cell 1): 0.416509
#> theta (batch 1): 0.220284
#> betas:    1.2407  -0.2570   0.4831   0.6402   0.6126   0.5554   0.6291  -0.0693   0.6241  -0.2704  -0.1496   0.5193
#> 
#> sigma: 0.159976
#> lambda (gene 1): 1.49437
#> epsilon (gene 1): 0.205246
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 684.656
#> delta (gene 1): 0.978548
#> phi (cell 1): 0.811522
#> s (cell 1): 0.490377
#> nu (cell 1): 0.42519
#> theta (batch 1): 0.312248
#> betas:    0.8523  -0.0557   0.5465   0.9914   0.7109   0.0694   0.5916  -0.1101  -0.1657  -0.1365  -0.4386  -0.1753
#> 
#> sigma: 0.151128
#> lambda (gene 1): 0.950802
#> epsilon (gene 1): 0.128302
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 981.468
#> delta (gene 1): 0.660412
#> phi (cell 1): 0.862766
#> s (cell 1): 0.422613
#> nu (cell 1): 0.423217
#> theta (batch 1): 0.0122239
#> betas:    0.8538   0.0547   0.0837   0.5486   0.7984   0.6088  -0.2137  -0.0408  -0.0300  -0.8398  -0.5636  -0.7784
#> 
#> sigma: 0.194166
#> lambda (gene 1): 0.527181
#> epsilon (gene 1): -0.106506
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 856.213
#> delta (gene 1): 1.01921
#> phi (cell 1): 0.691317
#> s (cell 1): 0.40065
#> nu (cell 1): 0.434845
#> theta (batch 1): 0.0125897
#> betas:    0.9062  -0.0411   0.8795   0.4499   0.3424   0.6134   0.3052  -0.1147  -0.0872  -0.2959  -0.2122  -0.6612
#> 
#> sigma: 0.161443
#> lambda (gene 1): 1.02571
#> epsilon (gene 1): 0.104051
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 619.709
#> delta (gene 1): 0.804437
#> phi (cell 1): 0.75578
#> s (cell 1): 0.390004
#> nu (cell 1): 0.416778
#> theta (batch 1): 0.00300712
#> betas:    1.2799  -0.0979   0.3101   0.9050   0.7315  -0.0508   0.5894  -0.1019   0.0128  -0.5111  -0.0887  -0.7375
#> 
#> sigma: 0.141411
#> lambda (gene 1): 0.881682
#> epsilon (gene 1): -0.0617403
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 778.815
#> delta (gene 1): 0.575211
#> phi (cell 1): 0.743637
#> s (cell 1): 0.455719
#> nu (cell 1): 0.412926
#> theta (batch 1): 0.010494
#> betas:    0.9876  -0.0774   0.7924   1.2015   0.1917   0.5508   0.2120   0.0518  -0.0510  -0.4481  -0.1764  -0.5021
#> 
#> sigma: 0.158671
#> lambda (gene 1): 1.45089
#> epsilon (gene 1): -0.285478
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 656.659
#> delta (gene 1): 0.664049
#> phi (cell 1): 0.713279
#> s (cell 1): 0.459788
#> nu (cell 1): 0.416757
#> theta (batch 1): 0.0210308
#> betas:    0.9527  -0.0861   0.9753   0.9897  -0.4279   0.9956   0.2120   0.1648   0.0467   0.0332  -0.9106   0.3857
#> 
#> sigma: 0.148021
#> lambda (gene 1): 0.375604
#> epsilon (gene 1): -0.315319
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 652.603
#> delta (gene 1): 0.665955
#> phi (cell 1): 0.75656
#> s (cell 1): 0.622682
#> nu (cell 1): 0.440219
#> theta (batch 1): 0.0670853
#> betas:    0.4842  -0.0423   1.1572   0.2631   1.0726   0.3357   0.5208   0.2355  -0.1759  -0.1669  -0.2575  -0.1556
#> 
#> sigma: 0.159788
#> lambda (gene 1): 0.733498
#> epsilon (gene 1): -0.156013
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 723.476
#> delta (gene 1): 0.803566
#> phi (cell 1): 0.793329
#> s (cell 1): 0.873548
#> nu (cell 1): 0.415808
#> theta (batch 1): 0.186033
#> betas:    0.9893  -0.0451   0.4075   0.5839   0.7071   0.2143   0.4907  -0.4983   0.3261  -0.5318  -0.7499  -0.2058
#> 
#> sigma: 0.152074
#> lambda (gene 1): 2.12437
#> epsilon (gene 1): 0.129783
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 849.964
#> delta (gene 1): 0.681812
#> phi (cell 1): 0.689664
#> s (cell 1): 1.09233
#> nu (cell 1): 0.435632
#> theta (batch 1): 0.310191
#> betas:    0.4785  -0.0326   0.5299   1.0791   0.7228   0.9770  -0.5126   0.9763  -0.2128  -0.4697   0.0446  -0.3416
#> 
#> sigma: 0.201573
#> lambda (gene 1): 1.12493
#> epsilon (gene 1): -0.16482
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 615.531
#> delta (gene 1): 0.728051
#> phi (cell 1): 0.727999
#> s (cell 1): 0.363214
#> nu (cell 1): 0.426431
#> theta (batch 1): 0.117728
#> betas:    0.8484  -0.0294   0.3709   0.9342   0.3806   0.4126   0.7905  -0.3077  -0.1388  -0.0777  -0.6338  -0.4122
#> 
#> sigma: 0.195984
#> lambda (gene 1): 1.26372
#> epsilon (gene 1): -0.162389
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 604.291
#> delta (gene 1): 0.619398
#> phi (cell 1): 0.771269
#> s (cell 1): 0.237755
#> nu (cell 1): 0.427616
#> theta (batch 1): 0.204964
#> betas:    1.4960  -0.1701   0.9921   0.6543   0.5660   0.2406   0.4454  -0.2147   0.2368  -0.4352  -0.5777   0.2898
#> 
#> sigma: 0.184616
#> lambda (gene 1): 1.51069
#> epsilon (gene 1): -0.178465
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 666.26
#> delta (gene 1): 0.668427
#> phi (cell 1): 0.712249
#> s (cell 1): 0.574482
#> nu (cell 1): 0.429321
#> theta (batch 1): 0.385703
#> betas:    1.2097  -0.2017   0.8234   0.9452   0.5648   0.3109   0.7425   0.0673  -0.1734   0.2178  -0.2242  -0.0387
#> 
#> sigma: 0.194398
#> lambda (gene 1): 0.371528
#> epsilon (gene 1): -0.233061
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 614.665
#> delta (gene 1): 0.521085
#> phi (cell 1): 0.699751
#> s (cell 1): 0.58705
#> nu (cell 1): 0.426601
#> theta (batch 1): 0.3711
#> betas:    0.7705  -0.0165   0.5073   0.9079   0.7190   0.0648   0.5920   0.2784  -0.6166   0.0793  -0.7192  -0.4674
#> 
#> sigma: 0.116885
#> lambda (gene 1): 2.20013
#> epsilon (gene 1): -0.385982
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 644.217
#> delta (gene 1): 0.616633
#> phi (cell 1): 0.732593
#> s (cell 1): 0.904392
#> nu (cell 1): 0.411641
#> theta (batch 1): 0.267206
#> betas:    0.8095  -0.0190   0.3352   1.2458   0.5563   0.4701   0.5727  -0.5525   0.3828  -0.4669  -0.7255  -0.3131
#> 
#> sigma: 0.144248
#> lambda (gene 1): 2.70201
#> epsilon (gene 1): -0.210027
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 822.702
#> delta (gene 1): 0.621613
#> phi (cell 1): 0.765394
#> s (cell 1): 0.712158
#> nu (cell 1): 0.429369
#> theta (batch 1): 0.281712
#> betas:    0.5935  -0.0760   0.6720   0.5528   0.4629   1.0300   0.5044  -0.3688   0.6254   0.0108  -0.8653   0.1694
#> 
#> sigma: 0.116379
#> lambda (gene 1): 1.14618
#> epsilon (gene 1): -0.311774
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 779.077
#> delta (gene 1): 0.848817
#> phi (cell 1): 0.720222
#> s (cell 1): 0.389959
#> nu (cell 1): 0.42311
#> theta (batch 1): 0.187951
#> betas:    1.2230  -0.0814   0.4177   0.4219   0.7083   0.4081   0.0789   0.2134  -0.4514  -0.3255  -0.2975  -0.3024
#> 
#> sigma: 0.176483
#> lambda (gene 1): 2.91321
#> epsilon (gene 1): 0.037537
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 830.088
#> delta (gene 1): 0.859725
#> phi (cell 1): 0.711953
#> s (cell 1): 1.00637
#> nu (cell 1): 0.435956
#> theta (batch 1): 0.247433
#> betas:    0.8637  -0.1432   0.9070   0.1500   1.0921   0.6234   0.3432   0.1381   0.1805   0.1394  -0.6779   0.2316
#> 
#> sigma: 0.140731
#> lambda (gene 1): 0.826931
#> epsilon (gene 1): 0.0887455
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 690.868
#> delta (gene 1): 0.935902
#> phi (cell 1): 0.65033
#> s (cell 1): 0.722295
#> nu (cell 1): 0.423406
#> theta (batch 1): 0.376187
#> betas:    1.2262  -0.0841   0.4401   0.3892   0.5692   0.6687   0.2341  -0.1556  -0.0683  -0.3626  -0.4277  -0.4311
#> 
#> sigma: 0.109784
#> lambda (gene 1): 1.69336
#> epsilon (gene 1): 0.147166
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 667.998
#> delta (gene 1): 0.709578
#> phi (cell 1): 0.730601
#> s (cell 1): 0.842142
#> nu (cell 1): 0.422413
#> theta (batch 1): 0.144403
#> betas:    1.1579  -0.1282   0.6812   0.4833   0.6820   0.5179  -0.1191   0.5929  -0.4524   0.0793  -0.3848  -0.5289
#> 
#> sigma: 0.186683
#> lambda (gene 1): 0.303599
#> epsilon (gene 1): -0.1082
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 668.99
#> delta (gene 1): 0.593094
#> phi (cell 1): 0.653739
#> s (cell 1): 0.435488
#> nu (cell 1): 0.415217
#> theta (batch 1): 0.141933
#> betas:    0.4873  -0.0635   1.5765   0.2485   1.1271   0.7992  -0.0763   0.6396  -0.1549   0.0140  -0.3141  -0.0949
#> 
#> sigma: 0.161833
#> lambda (gene 1): 0.553376
#> epsilon (gene 1): -0.368284
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 763.062
#> delta (gene 1): 0.643809
#> phi (cell 1): 0.707309
#> s (cell 1): 1.05886
#> nu (cell 1): 0.405776
#> theta (batch 1): 0.231292
#> betas:    1.1195  -0.1290   1.2394   0.6532  -0.0664   0.8623   0.2488  -0.0017   0.1277  -0.3754   0.0357  -0.6906
#> 
#> sigma: 0.136749
#> lambda (gene 1): 1.19926
#> epsilon (gene 1): -0.247131
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 666.879
#> delta (gene 1): 0.539943
#> phi (cell 1): 0.717662
#> s (cell 1): 0.536781
#> nu (cell 1): 0.420319
#> theta (batch 1): 0.0494653
#> betas:    0.8704  -0.1309   0.6383   0.4486   0.8497   0.4747   0.2806   0.1847   0.0891  -0.1729  -0.0560  -0.3162
#> 
#> sigma: 0.132282
#> lambda (gene 1): 1.01417
#> epsilon (gene 1): -0.416713
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 668.865
#> delta (gene 1): 0.834737
#> phi (cell 1): 0.647504
#> s (cell 1): 0.396402
#> nu (cell 1): 0.418941
#> theta (batch 1): 0.0700117
#> betas:    1.3776  -0.0800   0.8736   0.6142  -0.0908   0.9009   0.0868  -0.2659  -0.0625  -0.5111  -0.6666  -0.2102
#> 
#> sigma: 0.183011
#> lambda (gene 1): 0.831292
#> epsilon (gene 1): 0.143387
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 807.644
#> delta (gene 1): 0.706177
#> phi (cell 1): 0.710551
#> s (cell 1): 0.60955
#> nu (cell 1): 0.424119
#> theta (batch 1): 0.0553022
#> betas:    0.5295  -0.0003   0.4127   1.4432   0.2219   0.5823   0.1127   0.3271   0.1281  -0.6801  -0.5209  -0.0311
#> 
#> sigma: 0.169263
#> lambda (gene 1): 1.93498
#> epsilon (gene 1): 0.0501788
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 572.395
#> delta (gene 1): 0.632648
#> phi (cell 1): 0.715162
#> s (cell 1): 0.425128
#> nu (cell 1): 0.427026
#> theta (batch 1): 0.157347
#> betas:    0.8674  -0.1436   0.7785   0.9374   0.6255   0.1266   0.7723   0.4981  -0.0915  -0.2899  -0.0486   0.0036
#> 
#> sigma: 0.117005
#> lambda (gene 1): 0.78466
#> epsilon (gene 1): -0.118148
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 597.595
#> delta (gene 1): 0.463964
#> phi (cell 1): 0.663353
#> s (cell 1): 1.61622
#> nu (cell 1): 0.418714
#> theta (batch 1): 0.189191
#> betas:    0.7737  -0.1284   0.6421   1.2309   0.5024   0.3833   0.6336   0.0730   0.1922  -0.3580  -0.0917  -0.0650
#> 
#> sigma: 0.164267
#> lambda (gene 1): 0.0692452
#> epsilon (gene 1): -0.405325
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 711.54
#> delta (gene 1): 0.670173
#> phi (cell 1): 0.655755
#> s (cell 1): 1.31482
#> nu (cell 1): 0.422713
#> theta (batch 1): 0.764317
#> betas:    1.2412  -0.1886   0.5387   0.8630   0.3004   0.4836   0.6303  -0.2547   0.5921  -0.6681   0.0360  -0.1082
#> 
#> sigma: 0.141832
#> lambda (gene 1): 1.85425
#> epsilon (gene 1): -0.0706517
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 622.498
#> delta (gene 1): 0.471971
#> phi (cell 1): 0.650427
#> s (cell 1): 0.29228
#> nu (cell 1): 0.43126
#> theta (batch 1): 0.221446
#> betas:    0.7665  -0.1223   0.5051   0.9441   0.6529   0.2129   0.9109   0.1163   0.2969  -0.6142  -0.3563   0.5968
#> 
#> sigma: 0.203345
#> lambda (gene 1): 1.97518
#> epsilon (gene 1): -0.243721
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 737.943
#> delta (gene 1): 0.729999
#> phi (cell 1): 0.730351
#> s (cell 1): 0.410059
#> nu (cell 1): 0.422474
#> theta (batch 1): 0.163059
#> betas:    0.2500   0.0051   1.1841   0.5425   0.8636   1.0525   0.0055   0.2301   0.2146  -0.4677  -0.0929  -0.4596
#> 
#> sigma: 0.178472
#> lambda (gene 1): 1.42771
#> epsilon (gene 1): -0.11915
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 658.276
#> delta (gene 1): 0.689068
#> phi (cell 1): 0.748042
#> s (cell 1): 0.40056
#> nu (cell 1): 0.414179
#> theta (batch 1): 0.399834
#> betas:    1.1352  -0.0888   0.4688   0.8362   0.5855   0.4336   0.0826  -0.0238   0.0548  -0.3276  -0.3904  -0.5091
#> 
#> sigma: 0.138406
#> lambda (gene 1): 1.73694
#> epsilon (gene 1): -0.179254
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 836.904
#> delta (gene 1): 0.716496
#> phi (cell 1): 0.692299
#> s (cell 1): 0.561857
#> nu (cell 1): 0.433085
#> theta (batch 1): 0.273044
#> betas:    1.0818  -0.0989   0.8898   0.3655   0.7119   0.4581   0.1539   0.2294  -0.1627  -0.2700  -0.4085  -0.2170
#> 
#> sigma: 0.138179
#> lambda (gene 1): 1.60066
#> epsilon (gene 1): -0.0628661
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 684.119
#> delta (gene 1): 0.441274
#> phi (cell 1): 0.715171
#> s (cell 1): 0.429972
#> nu (cell 1): 0.42074
#> theta (batch 1): 0.378796
#> betas:    0.5105  -0.0488   1.2751   0.3383   0.8258   0.7039   0.6210   0.0050  -0.1224   0.1074  -0.6530  -0.2950
#> 
#> sigma: 0.205127
#> lambda (gene 1): 0.185945
#> epsilon (gene 1): -0.472015
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 692.805
#> delta (gene 1): 0.912277
#> phi (cell 1): 0.660895
#> s (cell 1): 0.544956
#> nu (cell 1): 0.435159
#> theta (batch 1): 0.131464
#> betas:    0.5176   0.0366   0.6708   0.8233   0.7586   0.4666  -0.0681   0.2516   0.0410  -0.4542  -0.6741  -0.4213
#> 
#> sigma: 0.139179
#> lambda (gene 1): 0.55585
#> epsilon (gene 1): 0.157103
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 760.681
#> delta (gene 1): 0.658306
#> phi (cell 1): 0.701036
#> s (cell 1): 0.550356
#> nu (cell 1): 0.427186
#> theta (batch 1): 0.0832256
#> betas:    1.2807  -0.1923   0.6859   0.7306   0.3512   0.5087   0.6415  -0.2038   0.1914  -0.1991  -0.0735  -0.2876
#> 
#> sigma: 0.190468
#> lambda (gene 1): 0.841005
#> epsilon (gene 1): -0.185783
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 638.333
#> delta (gene 1): 0.893873
#> phi (cell 1): 0.634948
#> s (cell 1): 0.690865
#> nu (cell 1): 0.413727
#> theta (batch 1): 0.136353
#> betas:    0.4973  -0.0686   0.2318   1.0820   0.7130   0.7353   0.3433   0.1130   0.0897   0.1488  -0.6558   0.1007
#> 
#> sigma: 0.178659
#> lambda (gene 1): 0.830847
#> epsilon (gene 1): 0.073822
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 645.973
#> delta (gene 1): 0.814532
#> phi (cell 1): 0.693124
#> s (cell 1): 0.237949
#> nu (cell 1): 0.436831
#> theta (batch 1): 0.167058
#> betas:    0.5373   0.0254   0.6932   1.0884   0.5677   0.4576   0.3522   0.0864  -0.1987  -0.2917  -0.8096  -0.0106
#> 
#> sigma: 0.131488
#> lambda (gene 1): 0.807971
#> epsilon (gene 1): 0.0966618
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 553.89
#> delta (gene 1): 0.715188
#> phi (cell 1): 0.714831
#> s (cell 1): 0.514264
#> nu (cell 1): 0.417924
#> theta (batch 1): 0.0464578
#> betas:    1.2355  -0.2645   0.8825   0.3972   0.7045   0.1342   1.1684  -0.1082   0.0995   0.0494   0.2541   0.1387
#> 
#> sigma: 0.159969
#> lambda (gene 1): 0.995701
#> epsilon (gene 1): -0.176508
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 641.843
#> delta (gene 1): 0.504669
#> phi (cell 1): 0.64875
#> s (cell 1): 0.426145
#> nu (cell 1): 0.430297
#> theta (batch 1): 0.186613
#> betas:    1.1821  -0.1383   0.3356   0.4718   0.7589   0.3893   0.1140   0.0903  -0.0205  -0.0121  -0.5091   0.0369
#> 
#> sigma: 0.165967
#> lambda (gene 1): 0.206109
#> epsilon (gene 1): -0.567145
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 748.45
#> delta (gene 1): 0.900544
#> phi (cell 1): 0.689884
#> s (cell 1): 0.862938
#> nu (cell 1): 0.426112
#> theta (batch 1): 0.335555
#> betas:    0.4850  -0.0127   1.1010   0.5352   1.1610   0.1449   0.1904   0.6487  -0.2882  -0.2038  -0.6528   0.0117
#> 
#> sigma: 0.139823
#> lambda (gene 1): 1.5005
#> epsilon (gene 1): 0.235252
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 609.922
#> delta (gene 1): 0.594388
#> phi (cell 1): 0.670334
#> s (cell 1): 0.355219
#> nu (cell 1): 0.421918
#> theta (batch 1): 0.179532
#> betas:    0.6876  -0.1984   1.2005   0.5649   1.0994   0.6785   0.5718   0.5038  -0.0200   0.1808   0.2198   0.2026
#> 
#> sigma: 0.159074
#> lambda (gene 1): 0.667209
#> epsilon (gene 1): -0.409474
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 745.491
#> delta (gene 1): 1.00398
#> phi (cell 1): 0.722907
#> s (cell 1): 0.612503
#> nu (cell 1): 0.421992
#> theta (batch 1): 0.216989
#> betas:    0.7132  -0.0789   0.6323   0.4663   1.0163   0.3295   0.4664   0.0228  -0.2161   0.2659  -0.6912   0.1299
#> 
#> sigma: 0.13544
#> lambda (gene 1): 1.72925
#> epsilon (gene 1): 0.158802
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 666.442
#> delta (gene 1): 0.510681
#> phi (cell 1): 0.68988
#> s (cell 1): 1.98756
#> nu (cell 1): 0.428771
#> theta (batch 1): 0.234197
#> betas:    0.3893  -0.0131   1.0035   0.6815   0.7951   0.6337   0.5253  -0.0815   0.4892  -0.3104  -0.9100   0.5868
#> 
#> sigma: 0.194759
#> lambda (gene 1): 0.962071
#> epsilon (gene 1): -0.462308
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 753.045
#> delta (gene 1): 0.705194
#> phi (cell 1): 0.644058
#> s (cell 1): 0.861383
#> nu (cell 1): 0.424977
#> theta (batch 1): 0.2794
#> betas:    0.9545  -0.0269   0.0287   1.2091   0.0023   0.7389   0.2935  -0.1468   0.1614  -0.4045  -0.7145  -0.3300
#> 
#> sigma: 0.156205
#> lambda (gene 1): 1.99828
#> epsilon (gene 1): -0.166803
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 746.484
#> delta (gene 1): 0.603808
#> phi (cell 1): 0.662929
#> s (cell 1): 1.4037
#> nu (cell 1): 0.421217
#> theta (batch 1): 0.424729
#> betas:    0.7464  -0.0722   1.0074   0.4192   0.6383   0.4117   0.5609   0.0396   0.0854  -0.2184  -0.3980  -0.3210
#> 
#> sigma: 0.12683
#> lambda (gene 1): 0.478327
#> epsilon (gene 1): -0.240226
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 746.421
#> delta (gene 1): 0.721223
#> phi (cell 1): 0.676942
#> s (cell 1): 0.529007
#> nu (cell 1): 0.412012
#> theta (batch 1): 0.237775
#> betas:    0.5988  -0.0369   1.3729   0.4077   0.7579   0.5909   0.1195   0.4286  -0.1902  -0.3988  -0.3169  -0.3454
#> 
#> sigma: 0.15326
#> lambda (gene 1): 1.49284
#> epsilon (gene 1): 0.0862169
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 808.449
#> delta (gene 1): 0.502027
#> phi (cell 1): 0.730196
#> s (cell 1): 0.539999
#> nu (cell 1): 0.423761
#> theta (batch 1): 0.242813
#> betas:    0.8129   0.0653   0.5551   0.7691   0.3397   0.9114  -0.1498  -0.0425  -0.2669  -0.5361  -0.9276  -0.7877
#> 
#> sigma: 0.17992
#> lambda (gene 1): 0.549877
#> epsilon (gene 1): -0.28112
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 699.207
#> delta (gene 1): 0.708572
#> phi (cell 1): 0.643928
#> s (cell 1): 1.2425
#> nu (cell 1): 0.403668
#> theta (batch 1): 0.584735
#> betas:    0.8590   0.0559   1.1200   0.8472  -0.1274   0.7368   0.2086  -0.4915   0.1740  -0.7028  -0.4606  -1.0154
#> 
#> sigma: 0.224414
#> lambda (gene 1): 1.27753
#> epsilon (gene 1): -0.246907
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 7980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 595.062
#> delta (gene 1): 0.845787
#> phi (cell 1): 0.658392
#> s (cell 1): 0.604752
#> nu (cell 1): 0.437569
#> theta (batch 1): 0.246851
#> betas:    1.0413  -0.0943   0.7209   0.6720   0.7047   0.5453   0.1107  -0.0715   0.3953  -0.4621  -0.3118  -0.2402
#> 
#> sigma: 0.188368
#> lambda (gene 1): 2.90339
#> epsilon (gene 1): -0.0601985
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 513.696
#> delta (gene 1): 0.58776
#> phi (cell 1): 0.717253
#> s (cell 1): 0.328537
#> nu (cell 1): 0.425486
#> theta (batch 1): 0.116087
#> betas:    0.8199  -0.0816   0.0190   0.3496   1.3431   0.2883   0.4034   0.1095   0.0353  -0.1516  -0.7620   0.1773
#> 
#> sigma: 0.162838
#> lambda (gene 1): 1.23373
#> epsilon (gene 1): -0.189134
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 862.24
#> delta (gene 1): 0.698942
#> phi (cell 1): 0.809396
#> s (cell 1): 0.499746
#> nu (cell 1): 0.422248
#> theta (batch 1): 0.0877749
#> betas:    1.3662  -0.0937   0.4777   0.4389   0.3602   0.1730   0.5032  -0.1384   0.0871  -0.8962  -0.1643  -0.4112
#> 
#> sigma: 0.174173
#> lambda (gene 1): 0.554288
#> epsilon (gene 1): -0.00896677
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 592.18
#> delta (gene 1): 0.843854
#> phi (cell 1): 0.78388
#> s (cell 1): 0.594463
#> nu (cell 1): 0.411974
#> theta (batch 1): 0.048502
#> betas:    0.7153  -0.1427   0.7127   1.2431   0.1832   0.9878   0.2757   0.2919   0.2898  -0.2124  -0.3211   0.4472
#> 
#> sigma: 0.143343
#> lambda (gene 1): 0.680102
#> epsilon (gene 1): 0.128193
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 671.171
#> delta (gene 1): 0.76713
#> phi (cell 1): 0.770213
#> s (cell 1): 0.543615
#> nu (cell 1): 0.413407
#> theta (batch 1): 0.0182874
#> betas:    0.4279   0.0022   0.6689   1.3170   0.5862   0.5542   0.5018  -0.4599   0.3972  -0.2904  -0.4901  -0.3995
#> 
#> sigma: 0.160261
#> lambda (gene 1): 0.615588
#> epsilon (gene 1): -0.107975
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 805.65
#> delta (gene 1): 0.708989
#> phi (cell 1): 0.758449
#> s (cell 1): 0.45549
#> nu (cell 1): 0.426727
#> theta (batch 1): 0.0120523
#> betas:    0.8455  -0.1170   0.2799   1.2899   0.6068   0.3604   0.4038   0.2934  -0.1200   0.1448  -0.3965  -0.2841
#> 
#> sigma: 0.15911
#> lambda (gene 1): 0.257447
#> epsilon (gene 1): -0.168623
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 543.339
#> delta (gene 1): 0.78058
#> phi (cell 1): 0.831059
#> s (cell 1): 0.357039
#> nu (cell 1): 0.404272
#> theta (batch 1): 0.00843239
#> betas:    1.2816  -0.3449   0.6906   0.4891   0.4292   0.9508   0.3117   0.4767   0.4312   0.2139  -0.0849   0.7498
#> 
#> sigma: 0.246117
#> lambda (gene 1): 1.10819
#> epsilon (gene 1): 0.00531342
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 569.974
#> delta (gene 1): 0.595091
#> phi (cell 1): 0.786537
#> s (cell 1): 0.326556
#> nu (cell 1): 0.408464
#> theta (batch 1): 0.0395067
#> betas:    0.5384  -0.1036   0.8286   0.5720   0.7949   1.2426  -0.1475   0.4974   0.4762  -0.3961   0.0115   0.1513
#> 
#> sigma: 0.182659
#> lambda (gene 1): 0.581209
#> epsilon (gene 1): -0.43936
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 489.998
#> delta (gene 1): 0.507742
#> phi (cell 1): 0.734602
#> s (cell 1): 0.348255
#> nu (cell 1): 0.418523
#> theta (batch 1): 0.0901981
#> betas:    1.7603  -0.2000   0.6675   0.4454   0.3168   0.7174  -0.0183  -0.1089   0.0042  -0.0719  -0.5885  -0.2198
#> 
#> sigma: 0.189205
#> lambda (gene 1): 0.974773
#> epsilon (gene 1): -0.504033
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 735.52
#> delta (gene 1): 0.784457
#> phi (cell 1): 0.727775
#> s (cell 1): 0.463587
#> nu (cell 1): 0.420113
#> theta (batch 1): 0.0192497
#> betas:    0.6822   0.0308   0.5109   1.3149   0.1476   0.5643   0.3786  -0.3147   0.1149  -0.2650  -0.9581  -0.6724
#> 
#> sigma: 0.152518
#> lambda (gene 1): 2.98883
#> epsilon (gene 1): 0.0383106
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 880.937
#> delta (gene 1): 0.798565
#> phi (cell 1): 0.71104
#> s (cell 1): 0.402462
#> nu (cell 1): 0.427599
#> theta (batch 1): 0.0999073
#> betas:    1.1494  -0.1171   0.2202   0.7682   0.4975   0.6929   0.0046  -0.1452   0.2552  -0.2002  -0.5831  -0.1756
#> 
#> sigma: 0.145224
#> lambda (gene 1): 0.200675
#> epsilon (gene 1): -0.0371061
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 630.521
#> delta (gene 1): 1.08596
#> phi (cell 1): 0.681184
#> s (cell 1): 0.579939
#> nu (cell 1): 0.427663
#> theta (batch 1): 0.105184
#> betas:    1.0535  -0.0627   0.0665   0.5565   0.9043   0.2548   0.2803   0.2886  -0.3745  -0.2464  -0.3347  -0.6625
#> 
#> sigma: 0.151932
#> lambda (gene 1): 0.8932
#> epsilon (gene 1): 0.311751
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 632.556
#> delta (gene 1): 1.35824
#> phi (cell 1): 0.704017
#> s (cell 1): 0.501226
#> nu (cell 1): 0.406323
#> theta (batch 1): 0.202572
#> betas:    0.6560  -0.2179   0.8584   0.5151   1.1900   0.4304   0.7397   0.1614   0.6379  -0.3200   0.7138   0.1765
#> 
#> sigma: 0.177817
#> lambda (gene 1): 1.08056
#> epsilon (gene 1): 0.348463
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 672.313
#> delta (gene 1): 0.681549
#> phi (cell 1): 0.808552
#> s (cell 1): 0.303354
#> nu (cell 1): 0.409454
#> theta (batch 1): 0.18844
#> betas:    0.5614  -0.0003   0.7366   0.7495   0.8302   0.4379   0.1017   0.4861  -0.3036  -0.3994  -0.3139  -0.4440
#> 
#> sigma: 0.145026
#> lambda (gene 1): 0.35751
#> epsilon (gene 1): -0.0870837
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 648.346
#> delta (gene 1): 0.702525
#> phi (cell 1): 0.825395
#> s (cell 1): 0.812951
#> nu (cell 1): 0.413984
#> theta (batch 1): 0.125541
#> betas:    1.0924  -0.0911   0.6779   0.4663   0.8321   0.2834   0.2752  -0.3034   0.4024  -0.5494  -0.6036   0.2400
#> 
#> sigma: 0.173402
#> lambda (gene 1): 1.76458
#> epsilon (gene 1): -0.109782
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 933.2
#> delta (gene 1): 0.742343
#> phi (cell 1): 0.747485
#> s (cell 1): 0.810098
#> nu (cell 1): 0.424872
#> theta (batch 1): 0.38606
#> betas:    0.8528  -0.0535   0.7910   0.6495   0.4321   0.8895  -0.1852   0.1500   0.0062   0.1744  -0.7069  -0.8192
#> 
#> sigma: 0.212346
#> lambda (gene 1): 2.5875
#> epsilon (gene 1): -0.221993
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 656.951
#> delta (gene 1): 0.479939
#> phi (cell 1): 0.701553
#> s (cell 1): 0.575484
#> nu (cell 1): 0.422892
#> theta (batch 1): 0.392941
#> betas:    0.6177  -0.0183   0.7141   0.6026   0.0972   1.0949   0.6013  -0.8465   0.8405  -0.4980  -0.5624   0.0312
#> 
#> sigma: 0.174813
#> lambda (gene 1): 0.560279
#> epsilon (gene 1): -0.704555
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 585.712
#> delta (gene 1): 0.699617
#> phi (cell 1): 0.74394
#> s (cell 1): 0.366273
#> nu (cell 1): 0.418597
#> theta (batch 1): 0.292723
#> betas:    1.0220  -0.0631   0.3489   0.8985   0.6668   0.5493  -0.2963   0.3733  -0.1264  -0.5767  -0.4853   0.1459
#> 
#> sigma: 0.162488
#> lambda (gene 1): 1.28396
#> epsilon (gene 1): -0.0459709
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 690.349
#> delta (gene 1): 0.887549
#> phi (cell 1): 0.737995
#> s (cell 1): 0.224047
#> nu (cell 1): 0.419422
#> theta (batch 1): 0.379469
#> betas:    0.5341   0.0547   0.9808   0.9529   0.2289   0.2508   0.6737  -0.2384  -0.2391   0.1685  -1.2962  -0.3910
#> 
#> sigma: 0.107741
#> lambda (gene 1): 0.942074
#> epsilon (gene 1): 0.0993326
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 959.625
#> delta (gene 1): 0.804974
#> phi (cell 1): 0.75037
#> s (cell 1): 0.502867
#> nu (cell 1): 0.415157
#> theta (batch 1): 0.262941
#> betas:    0.8796  -0.1014   0.9774   0.7156   0.3597   0.7828   0.1710   0.1745  -0.1387   0.0985  -0.6637   0.1539
#> 
#> sigma: 0.181592
#> lambda (gene 1): 1.51412
#> epsilon (gene 1): -0.0136238
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 814.129
#> delta (gene 1): 0.771077
#> phi (cell 1): 0.769654
#> s (cell 1): 0.339734
#> nu (cell 1): 0.419714
#> theta (batch 1): 0.154792
#> betas:    0.5960  -0.0650   0.5709   1.4584   0.3554   0.9923  -0.0136   0.3327  -0.1260   0.0701  -0.4675  -0.1814
#> 
#> sigma: 0.17572
#> lambda (gene 1): 1.03292
#> epsilon (gene 1): -0.07456
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 617.218
#> delta (gene 1): 0.611817
#> phi (cell 1): 0.72524
#> s (cell 1): 0.29149
#> nu (cell 1): 0.427346
#> theta (batch 1): 0.135121
#> betas:    0.5841  -0.0717   0.7019   0.7868   0.9999   0.6961  -0.0727   0.6346   0.0541  -0.5308  -0.1767   0.0754
#> 
#> sigma: 0.123123
#> lambda (gene 1): 0.474629
#> epsilon (gene 1): -0.125591
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1015.57
#> delta (gene 1): 0.562025
#> phi (cell 1): 0.794301
#> s (cell 1): 1.54666
#> nu (cell 1): 0.41891
#> theta (batch 1): 0.317898
#> betas:    1.0852  -0.1862   0.7459   0.8275   0.3435   0.6822   0.4330   0.0835   0.4534  -0.2180  -0.4555   0.5175
#> 
#> sigma: 0.15137
#> lambda (gene 1): 0.541588
#> epsilon (gene 1): -0.337382
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 818.735
#> delta (gene 1): 0.754628
#> phi (cell 1): 0.704542
#> s (cell 1): 0.340096
#> nu (cell 1): 0.409821
#> theta (batch 1): 0.401853
#> betas:    0.7626  -0.1939   0.4952   0.6917   0.9105   0.3271   0.8499   0.1586   0.3719   0.3615  -0.3775   0.1995
#> 
#> sigma: 0.228623
#> lambda (gene 1): 0.815189
#> epsilon (gene 1): -0.189191
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 603.606
#> delta (gene 1): 0.625864
#> phi (cell 1): 0.786023
#> s (cell 1): 0.256477
#> nu (cell 1): 0.411194
#> theta (batch 1): 0.371645
#> betas:    0.8911  -0.0554   0.9611   0.6166   0.7458   0.0306   0.3542   0.1911  -0.0288  -0.4196  -0.3417   0.0375
#> 
#> sigma: 0.133347
#> lambda (gene 1): 1.59874
#> epsilon (gene 1): -0.357516
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 658.416
#> delta (gene 1): 0.708833
#> phi (cell 1): 0.806394
#> s (cell 1): 0.29757
#> nu (cell 1): 0.431215
#> theta (batch 1): 0.330722
#> betas:    0.9051  -0.1492   0.6397   0.3905   0.8040   0.3414   0.7073   0.1077   0.1215  -0.2230   0.2231  -0.5693
#> 
#> sigma: 0.1523
#> lambda (gene 1): 1.9624
#> epsilon (gene 1): -0.165154
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 752.342
#> delta (gene 1): 0.768548
#> phi (cell 1): 0.733668
#> s (cell 1): 0.844894
#> nu (cell 1): 0.423001
#> theta (batch 1): 0.559416
#> betas:    0.8773  -0.1481   0.5352   0.8146   0.4245   0.7304   0.4720  -0.1628   0.3921   0.0434  -0.5773   0.2004
#> 
#> sigma: 0.143369
#> lambda (gene 1): 1.16099
#> epsilon (gene 1): -0.0548073
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 550.816
#> delta (gene 1): 0.719284
#> phi (cell 1): 0.863236
#> s (cell 1): 0.452619
#> nu (cell 1): 0.417171
#> theta (batch 1): 0.0910025
#> betas:    0.7729  -0.0268   0.2609   0.8687   0.7477   0.6273   0.1587   0.2961  -0.5644  -0.0712  -0.6872  -0.1453
#> 
#> sigma: 0.182512
#> lambda (gene 1): 1.63522
#> epsilon (gene 1): 0.00975221
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 480.854
#> delta (gene 1): 0.965221
#> phi (cell 1): 0.933964
#> s (cell 1): 0.365212
#> nu (cell 1): 0.415946
#> theta (batch 1): 0.0380931
#> betas:    0.2876  -0.0769   0.9000   0.9074   0.8333   0.2581   0.8865  -0.0583   0.2892   0.1003  -0.5055   0.3624
#> 
#> sigma: 0.166217
#> lambda (gene 1): 1.28248
#> epsilon (gene 1): 0.162266
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 909.803
#> delta (gene 1): 0.525959
#> phi (cell 1): 0.976251
#> s (cell 1): 0.470952
#> nu (cell 1): 0.41313
#> theta (batch 1): 0.0163371
#> betas:    1.2058  -0.0802   0.3585   1.2134   0.3572  -0.1185   0.8278  -0.0850  -0.2143  -0.1954  -0.6814  -0.0982
#> 
#> sigma: 0.227253
#> lambda (gene 1): 1.65475
#> epsilon (gene 1): -0.442642
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 560.906
#> delta (gene 1): 0.593437
#> phi (cell 1): 0.969348
#> s (cell 1): 0.39145
#> nu (cell 1): 0.417264
#> theta (batch 1): 0.0280694
#> betas:    0.7647  -0.0410   0.4199   0.8181   0.4651   0.5970   0.0538   0.2984   0.2071  -0.7358  -0.4312   0.1622
#> 
#> sigma: 0.160259
#> lambda (gene 1): 0.189762
#> epsilon (gene 1): -0.200415
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 778.191
#> delta (gene 1): 0.489865
#> phi (cell 1): 0.914885
#> s (cell 1): 0.430538
#> nu (cell 1): 0.420396
#> theta (batch 1): 0.00453248
#> betas:    0.3753   0.0461   0.7446   1.1738   0.6858   0.2953   0.6597   0.1143  -0.0826  -0.0362  -0.8174  -0.3127
#> 
#> sigma: 0.204227
#> lambda (gene 1): 0.932638
#> epsilon (gene 1): -0.639065
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 738.795
#> delta (gene 1): 0.977335
#> phi (cell 1): 0.855014
#> s (cell 1): 0.3975
#> nu (cell 1): 0.429073
#> theta (batch 1): 0.00265508
#> betas:    1.4196  -0.2847   0.5215   0.3681   0.4396   0.6069   0.3884   0.4400  -0.1281  -0.1184   0.2074   0.3576
#> 
#> sigma: 0.14186
#> lambda (gene 1): 0.44921
#> epsilon (gene 1): 0.28254
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 788.483
#> delta (gene 1): 0.656708
#> phi (cell 1): 0.853815
#> s (cell 1): 0.362504
#> nu (cell 1): 0.425412
#> theta (batch 1): 0.014317
#> betas:    0.8403  -0.2297   0.4825   0.9511   0.8696   0.6948   0.5395   0.2423   0.5702  -0.3788   0.4551   0.1879
#> 
#> sigma: 0.179571
#> lambda (gene 1): 0.675365
#> epsilon (gene 1): -0.198511
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 718.433
#> delta (gene 1): 0.588422
#> phi (cell 1): 0.814058
#> s (cell 1): 0.427986
#> nu (cell 1): 0.425989
#> theta (batch 1): 0.00577617
#> betas:    0.3766   0.0065   0.6287   0.6926   0.6346   0.6342   0.5844  -0.0290  -0.3827   0.5584  -0.9610  -0.1594
#> 
#> sigma: 0.208859
#> lambda (gene 1): 0.514487
#> epsilon (gene 1): -0.493836
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 686.241
#> delta (gene 1): 0.546768
#> phi (cell 1): 0.777835
#> s (cell 1): 0.463535
#> nu (cell 1): 0.405392
#> theta (batch 1): 0.0074702
#> betas:    0.1867   0.0356   0.9506   1.2362   0.6504   0.6433   0.9321  -0.5350   0.2588  -0.4639  -0.0522  -0.8317
#> 
#> sigma: 0.217921
#> lambda (gene 1): 1.30934
#> epsilon (gene 1): -0.371932
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 749.926
#> delta (gene 1): 1.09183
#> phi (cell 1): 0.764658
#> s (cell 1): 0.47985
#> nu (cell 1): 0.430571
#> theta (batch 1): 0.0207426
#> betas:    0.7691   0.0126   0.8181   0.3258   0.5934   0.5461  -0.1804   0.7128  -0.4960  -0.5081   0.0482  -1.0847
#> 
#> sigma: 0.161627
#> lambda (gene 1): 0.729162
#> epsilon (gene 1): 0.196282
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 835.599
#> delta (gene 1): 0.617142
#> phi (cell 1): 0.809642
#> s (cell 1): 0.862862
#> nu (cell 1): 0.431027
#> theta (batch 1): 0.225196
#> betas:    0.6856  -0.0808   0.2129   0.8356   0.8302   0.6665   0.2745   0.2493   0.0763  -0.1372  -0.4284   0.0559
#> 
#> sigma: 0.164349
#> lambda (gene 1): 1.41434
#> epsilon (gene 1): -0.297287
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 776.755
#> delta (gene 1): 0.483408
#> phi (cell 1): 0.827535
#> s (cell 1): 0.657715
#> nu (cell 1): 0.42425
#> theta (batch 1): 0.164994
#> betas:    0.7905   0.0090   0.9291   0.6546   0.4508   0.5408   0.1243  -0.2170   0.1923  -0.8423  -0.1136  -0.5017
#> 
#> sigma: 0.139227
#> lambda (gene 1): 0.802208
#> epsilon (gene 1): -0.590561
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 800.213
#> delta (gene 1): 0.536872
#> phi (cell 1): 0.862365
#> s (cell 1): 0.894664
#> nu (cell 1): 0.422242
#> theta (batch 1): 0.260385
#> betas:    1.2735  -0.1802   0.5811   0.8827   0.3708   0.5919   0.2111   0.0180   0.3934  -0.3388  -0.3361   0.1297
#> 
#> sigma: 0.121349
#> lambda (gene 1): 3.12623
#> epsilon (gene 1): -0.408523
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 683.384
#> delta (gene 1): 0.73342
#> phi (cell 1): 0.737193
#> s (cell 1): 0.572918
#> nu (cell 1): 0.419665
#> theta (batch 1): 0.298892
#> betas:    0.6894   0.0842   1.0218   0.1221   0.7599   0.5163   0.0827  -0.0250  -0.3571  -0.4712  -0.7086  -0.7732
#> 
#> sigma: 0.143732
#> lambda (gene 1): 0.696158
#> epsilon (gene 1): -0.0920337
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 646.501
#> delta (gene 1): 0.920138
#> phi (cell 1): 0.788878
#> s (cell 1): 0.218193
#> nu (cell 1): 0.422049
#> theta (batch 1): 0.0720386
#> betas:    1.0904  -0.2331   0.5125   0.8479   0.7947   0.3755   0.6400   0.0349   0.2458   0.4254  -0.3887   0.4506
#> 
#> sigma: 0.147211
#> lambda (gene 1): 0.908987
#> epsilon (gene 1): -0.0710796
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 601.729
#> delta (gene 1): 0.641873
#> phi (cell 1): 0.789854
#> s (cell 1): 0.75729
#> nu (cell 1): 0.433
#> theta (batch 1): 0.0716757
#> betas:    1.2101  -0.1002   0.3575   0.5849   0.9340   0.4301  -0.1032   0.3144  -0.1428  -0.5433  -0.1911  -0.5306
#> 
#> sigma: 0.144335
#> lambda (gene 1): 0.4733
#> epsilon (gene 1): -0.117002
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 583.928
#> delta (gene 1): 0.749393
#> phi (cell 1): 0.886362
#> s (cell 1): 0.292648
#> nu (cell 1): 0.410397
#> theta (batch 1): 0.0785553
#> betas:    1.0263  -0.0607   0.4102   0.5222   1.0619   0.1023  -0.1022   0.4188  -0.1299  -0.5707  -0.1288  -0.6778
#> 
#> sigma: 0.197746
#> lambda (gene 1): 1.29723
#> epsilon (gene 1): -0.0476254
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 698.744
#> delta (gene 1): 0.534856
#> phi (cell 1): 0.848361
#> s (cell 1): 0.478031
#> nu (cell 1): 0.436146
#> theta (batch 1): 0.584762
#> betas:    1.0623  -0.1497   0.5023   0.6779   0.5412   0.8024   0.1430   0.4545  -0.4443   0.1958  -0.2796  -0.2733
#> 
#> sigma: 0.115874
#> lambda (gene 1): 1.47344
#> epsilon (gene 1): -0.408811
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 696.401
#> delta (gene 1): 0.686143
#> phi (cell 1): 0.935392
#> s (cell 1): 3.51927
#> nu (cell 1): 0.43001
#> theta (batch 1): 0.344757
#> betas:    0.8066  -0.0146   0.7242   0.5317   0.8458  -0.0620   0.6343  -0.1425  -0.3874   0.2604  -1.1352  -0.3837
#> 
#> sigma: 0.174403
#> lambda (gene 1): 1.37055
#> epsilon (gene 1): -0.108159
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 706.12
#> delta (gene 1): 0.771714
#> phi (cell 1): 0.918493
#> s (cell 1): 0.402472
#> nu (cell 1): 0.415483
#> theta (batch 1): 0.365851
#> betas:    0.7390  -0.0571   0.2944   0.5652   1.4563   0.3557  -0.0806   0.4195  -0.0145  -0.2306  -0.3272  -0.3721
#> 
#> sigma: 0.166017
#> lambda (gene 1): 1.18772
#> epsilon (gene 1): -0.108325
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 830.921
#> delta (gene 1): 0.492995
#> phi (cell 1): 0.731488
#> s (cell 1): 0.731904
#> nu (cell 1): 0.415789
#> theta (batch 1): 0.451008
#> betas:    0.5779  -0.0811   1.1104   0.7432   0.2309   0.8301   0.2910   0.4248  -0.0857  -0.1293  -0.1724  -0.2500
#> 
#> sigma: 0.155195
#> lambda (gene 1): 1.3951
#> epsilon (gene 1): -0.458875
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 707.278
#> delta (gene 1): 0.631461
#> phi (cell 1): 0.661554
#> s (cell 1): 0.508969
#> nu (cell 1): 0.417477
#> theta (batch 1): 0.295741
#> betas:    0.6407  -0.0640   0.1303   0.7983   0.8555   0.1756   0.6129   0.3805  -0.5289   0.2430  -0.5671  -0.2812
#> 
#> sigma: 0.15327
#> lambda (gene 1): 0.242561
#> epsilon (gene 1): -0.160657
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 597
#> delta (gene 1): 0.574142
#> phi (cell 1): 0.648534
#> s (cell 1): 0.311483
#> nu (cell 1): 0.41113
#> theta (batch 1): 0.240741
#> betas:    0.4675  -0.0933   0.4740   0.9168   1.0476   0.5773   0.3661  -0.0494   0.9324  -0.7815   0.3116  -0.3901
#> 
#> sigma: 0.156708
#> lambda (gene 1): 0.145768
#> epsilon (gene 1): -0.355902
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 8980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 669.705
#> delta (gene 1): 0.627702
#> phi (cell 1): 0.677674
#> s (cell 1): 0.363021
#> nu (cell 1): 0.444611
#> theta (batch 1): 0.465827
#> betas:    0.8429   0.0043   0.7297   0.2112   0.9378   0.3988   0.1366  -0.3858   0.6410  -0.9048  -0.6209  -0.4686
#> 
#> sigma: 0.163311
#> lambda (gene 1): 1.13438
#> epsilon (gene 1): -0.1665
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 691.559
#> delta (gene 1): 0.704493
#> phi (cell 1): 0.68222
#> s (cell 1): 0.41363
#> nu (cell 1): 0.425797
#> theta (batch 1): 0.056605
#> betas:    0.7793  -0.0528   0.8263   0.8533   0.4314   0.8576  -0.1057   0.2845   0.1885  -0.5898  -0.1047  -0.4660
#> 
#> sigma: 0.145082
#> lambda (gene 1): 1.21047
#> epsilon (gene 1): -0.156352
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 608.523
#> delta (gene 1): 0.655068
#> phi (cell 1): 0.592745
#> s (cell 1): 0.559793
#> nu (cell 1): 0.418269
#> theta (batch 1): 0.296185
#> betas:    1.4058  -0.1255   0.8664   0.4690   0.6104   0.1108   0.1552   0.2317  -0.5182  -0.1404   0.0828  -0.6831
#> 
#> sigma: 0.124771
#> lambda (gene 1): 0.640602
#> epsilon (gene 1): -0.436173
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 638.842
#> delta (gene 1): 0.615631
#> phi (cell 1): 0.670787
#> s (cell 1): 0.556657
#> nu (cell 1): 0.423576
#> theta (batch 1): 0.102118
#> betas:    0.4755   0.0442   1.2457   0.4375   0.7465   0.7737  -0.1152   0.2579  -0.2201  -0.1631  -0.8607  -0.3029
#> 
#> sigma: 0.129628
#> lambda (gene 1): 0.825979
#> epsilon (gene 1): -0.257857
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 628.083
#> delta (gene 1): 0.721778
#> phi (cell 1): 0.689335
#> s (cell 1): 0.421065
#> nu (cell 1): 0.421966
#> theta (batch 1): 0.168062
#> betas:    1.3701  -0.2795   0.2652   0.6807   0.8566   0.2824   0.8980   0.3484  -0.2496   0.1638   0.0628   0.4330
#> 
#> sigma: 0.215468
#> lambda (gene 1): 1.58086
#> epsilon (gene 1): -0.116984
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 619.59
#> delta (gene 1): 0.631297
#> phi (cell 1): 0.806295
#> s (cell 1): 0.664882
#> nu (cell 1): 0.426394
#> theta (batch 1): 0.404733
#> betas:    0.7087  -0.0742   0.5848   0.3586   1.0412   0.3362   0.1324   0.5753  -0.1084  -0.2264  -0.5636   0.3676
#> 
#> sigma: 0.129174
#> lambda (gene 1): 1.25151
#> epsilon (gene 1): -0.183647
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 631.149
#> delta (gene 1): 0.692428
#> phi (cell 1): 0.775881
#> s (cell 1): 1.0504
#> nu (cell 1): 0.417436
#> theta (batch 1): 0.285994
#> betas:    0.8634  -0.1750   0.2036   0.7119   0.8993   1.0071  -0.2180   0.4799   0.2886  -0.3032   0.0547   0.0058
#> 
#> sigma: 0.160732
#> lambda (gene 1): 1.91399
#> epsilon (gene 1): -0.0919525
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 621.783
#> delta (gene 1): 0.483467
#> phi (cell 1): 0.744312
#> s (cell 1): 0.747983
#> nu (cell 1): 0.427204
#> theta (batch 1): 0.494006
#> betas:    0.9818  -0.0472   0.6339   0.8085   0.4821   0.3881   0.2156   0.0710  -0.0030  -0.5960  -0.4127  -0.1430
#> 
#> sigma: 0.188758
#> lambda (gene 1): 0.706117
#> epsilon (gene 1): -0.466276
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 576.499
#> delta (gene 1): 0.760052
#> phi (cell 1): 0.703357
#> s (cell 1): 0.438703
#> nu (cell 1): 0.412718
#> theta (batch 1): 0.197192
#> betas:    0.4928  -0.0569   0.6964   0.6593   0.9869   0.5981   0.3724  -0.0890   0.3542  -0.2395  -0.5873   0.3926
#> 
#> sigma: 0.123869
#> lambda (gene 1): 1.79274
#> epsilon (gene 1): -0.0421634
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 739.409
#> delta (gene 1): 0.621408
#> phi (cell 1): 0.663133
#> s (cell 1): 0.807723
#> nu (cell 1): 0.423383
#> theta (batch 1): 0.322661
#> betas:    1.1047  -0.1526   0.0250   1.1118   0.4363   0.5853   0.3975   0.2383  -0.1757  -0.2085   0.2733  -0.5036
#> 
#> sigma: 0.137147
#> lambda (gene 1): 0.126245
#> epsilon (gene 1): -0.370408
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 620.037
#> delta (gene 1): 0.684923
#> phi (cell 1): 0.686002
#> s (cell 1): 0.592032
#> nu (cell 1): 0.40721
#> theta (batch 1): 0.576365
#> betas:    0.1404   0.0260   0.5714   1.0494   0.5052   1.1421   0.3243  -0.2381   0.8035  -0.9163   0.1200  -0.5052
#> 
#> sigma: 0.17677
#> lambda (gene 1): 1.1914
#> epsilon (gene 1): -0.209099
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 542.971
#> delta (gene 1): 0.808744
#> phi (cell 1): 0.81374
#> s (cell 1): 0.454009
#> nu (cell 1): 0.424097
#> theta (batch 1): 0.202549
#> betas:    1.0648  -0.1306   0.1287   0.9842   0.7368   0.5141   0.1199   0.5287  -0.0885  -0.6694   0.3578  -0.2987
#> 
#> sigma: 0.146494
#> lambda (gene 1): 1.44462
#> epsilon (gene 1): 0.0129815
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 792.348
#> delta (gene 1): 0.678307
#> phi (cell 1): 0.717843
#> s (cell 1): 1.87129
#> nu (cell 1): 0.430122
#> theta (batch 1): 0.646604
#> betas:    0.5739  -0.0207   0.8932   0.9085   0.3019   0.6763   0.4267  -0.3205   0.4472  -0.6221  -0.5356  -0.2413
#> 
#> sigma: 0.187374
#> lambda (gene 1): 1.19386
#> epsilon (gene 1): 0.0345148
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 722.196
#> delta (gene 1): 0.422243
#> phi (cell 1): 0.687442
#> s (cell 1): 0.380567
#> nu (cell 1): 0.425513
#> theta (batch 1): 0.287624
#> betas:    0.1117  -0.0085   0.7991   1.2556   0.3530   1.3347  -0.0609   0.5368   0.3210  -0.5728  -0.2026   0.0505
#> 
#> sigma: 0.165793
#> lambda (gene 1): 0.714896
#> epsilon (gene 1): -0.530816
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 562.02
#> delta (gene 1): 0.63002
#> phi (cell 1): 0.773557
#> s (cell 1): 0.253034
#> nu (cell 1): 0.421615
#> theta (batch 1): 0.14371
#> betas:    0.7830  -0.0907   0.5181   0.6619   0.8104   0.4227   0.4843  -0.2330   0.5364  -0.6142  -0.1452  -0.1103
#> 
#> sigma: 0.119379
#> lambda (gene 1): 1.46889
#> epsilon (gene 1): -0.203936
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 858.785
#> delta (gene 1): 0.595067
#> phi (cell 1): 0.731594
#> s (cell 1): 1.00081
#> nu (cell 1): 0.429532
#> theta (batch 1): 0.446834
#> betas:    0.4948  -0.0163   0.8848   0.5317   0.8404   0.4254   0.5355   0.0815   0.0179  -0.1316  -0.5648  -0.5332
#> 
#> sigma: 0.193551
#> lambda (gene 1): 1.20322
#> epsilon (gene 1): -0.236753
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 552.566
#> delta (gene 1): 0.590398
#> phi (cell 1): 0.698637
#> s (cell 1): 0.382011
#> nu (cell 1): 0.410535
#> theta (batch 1): 0.164889
#> betas:    1.0146  -0.0933   1.0716   0.8172   0.2929   0.6854   0.6552  -0.4458   0.3765  -0.2830  -0.3926  -0.3454
#> 
#> sigma: 0.133894
#> lambda (gene 1): 0.529346
#> epsilon (gene 1): -0.413676
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 692.512
#> delta (gene 1): 0.50997
#> phi (cell 1): 0.76844
#> s (cell 1): 0.724926
#> nu (cell 1): 0.42768
#> theta (batch 1): 0.238255
#> betas:    1.1031  -0.0939   0.6314   0.3182   0.5789   0.3277   0.4895   0.0829  -0.2034   0.5909  -1.5105   0.4006
#> 
#> sigma: 0.241839
#> lambda (gene 1): 1.34098
#> epsilon (gene 1): -0.608263
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 780.49
#> delta (gene 1): 0.607332
#> phi (cell 1): 0.850542
#> s (cell 1): 0.666029
#> nu (cell 1): 0.429717
#> theta (batch 1): 0.382735
#> betas:    0.9859  -0.0347   0.5386   1.2302   0.0238   0.6477   0.1016   0.0333  -0.1506  -0.5751  -0.1676  -0.8316
#> 
#> sigma: 0.235562
#> lambda (gene 1): 0.398671
#> epsilon (gene 1): -0.225351
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 771.543
#> delta (gene 1): 0.576094
#> phi (cell 1): 0.754835
#> s (cell 1): 0.740374
#> nu (cell 1): 0.422845
#> theta (batch 1): 0.338008
#> betas:    0.7182  -0.1155   0.9287   0.4659   0.9050   0.5732   0.4055  -0.0478   0.6294  -0.3319  -0.1913  -0.0790
#> 
#> sigma: 0.166895
#> lambda (gene 1): 1.62643
#> epsilon (gene 1): -0.415464
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 660.401
#> delta (gene 1): 0.669122
#> phi (cell 1): 0.677557
#> s (cell 1): 0.54631
#> nu (cell 1): 0.41339
#> theta (batch 1): 0.141496
#> betas:    0.7117  -0.0722   0.5178   0.5486   0.9172   0.5937   0.1236   0.2263   0.3073  -0.3711  -0.1994  -0.2960
#> 
#> sigma: 0.139616
#> lambda (gene 1): 0.746828
#> epsilon (gene 1): -0.278177
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 916.955
#> delta (gene 1): 1.06691
#> phi (cell 1): 0.641301
#> s (cell 1): 0.460626
#> nu (cell 1): 0.419333
#> theta (batch 1): 0.16842
#> betas:    0.7956  -0.1591   0.9725   0.9848   0.9759   0.0821   0.5286   0.4122   0.1026  -0.1116   0.0806   0.0041
#> 
#> sigma: 0.159892
#> lambda (gene 1): 1.70543
#> epsilon (gene 1): 0.204159
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 565.426
#> delta (gene 1): 0.924498
#> phi (cell 1): 0.639919
#> s (cell 1): 0.543003
#> nu (cell 1): 0.412838
#> theta (batch 1): 0.39132
#> betas:    1.2014  -0.0538   0.3782   0.9485   0.1253   0.3009   0.3235  -0.3139   0.1302  -0.4226  -0.5520  -0.5476
#> 
#> sigma: 0.139288
#> lambda (gene 1): 1.20135
#> epsilon (gene 1): 0.049061
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 524.414
#> delta (gene 1): 0.709027
#> phi (cell 1): 0.694454
#> s (cell 1): 0.421084
#> nu (cell 1): 0.427279
#> theta (batch 1): 0.271163
#> betas:    0.6833  -0.1369   0.7556   1.0158   0.4137   0.9177   0.4293  -0.0752   0.6132  -0.1858  -0.7274   0.7143
#> 
#> sigma: 0.185723
#> lambda (gene 1): 0.963937
#> epsilon (gene 1): 0.0154254
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 687.133
#> delta (gene 1): 0.905232
#> phi (cell 1): 0.701386
#> s (cell 1): 0.286372
#> nu (cell 1): 0.415462
#> theta (batch 1): 0.173867
#> betas:    0.4424  -0.0533   0.0070   0.9538   0.8002   1.0619  -0.1514   0.1832   0.5908  -0.4001  -0.3270  -0.0975
#> 
#> sigma: 0.139008
#> lambda (gene 1): 2.21678
#> epsilon (gene 1): 0.0686847
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 812.653
#> delta (gene 1): 0.506899
#> phi (cell 1): 0.755155
#> s (cell 1): 0.502061
#> nu (cell 1): 0.42359
#> theta (batch 1): 0.197964
#> betas:    1.1876  -0.0993   0.3939   0.7815   0.4805   0.6140  -0.0279  -0.3596   0.8885  -0.9975  -0.5196   0.2848
#> 
#> sigma: 0.191615
#> lambda (gene 1): 0.734147
#> epsilon (gene 1): -0.401393
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 941.488
#> delta (gene 1): 0.739799
#> phi (cell 1): 0.760968
#> s (cell 1): 0.842214
#> nu (cell 1): 0.428387
#> theta (batch 1): 0.3511
#> betas:    1.0558  -0.1653   1.0054   0.6384   0.5770   0.6762  -0.0979   0.3575   0.4582  -0.4046  -0.2265   0.1513
#> 
#> sigma: 0.121006
#> lambda (gene 1): 0.664065
#> epsilon (gene 1): -0.0975619
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 786.573
#> delta (gene 1): 0.729187
#> phi (cell 1): 0.741567
#> s (cell 1): 0.526068
#> nu (cell 1): 0.434291
#> theta (batch 1): 0.114026
#> betas:    0.6147  -0.0469   0.3642   0.4060   1.0700   0.6946   0.0833   0.1396   0.0752  -0.1054  -0.7353  -0.0867
#> 
#> sigma: 0.147368
#> lambda (gene 1): 0.956126
#> epsilon (gene 1): -0.0104771
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 652.838
#> delta (gene 1): 0.618855
#> phi (cell 1): 0.754334
#> s (cell 1): 0.33193
#> nu (cell 1): 0.42334
#> theta (batch 1): 0.209747
#> betas:    0.8147  -0.1158   0.6344   0.6000   0.7402   0.7883  -0.0668   0.8028  -0.3145   0.0316  -0.2386  -0.2434
#> 
#> sigma: 0.145978
#> lambda (gene 1): 1.38719
#> epsilon (gene 1): -0.267823
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 817.049
#> delta (gene 1): 0.980659
#> phi (cell 1): 0.754808
#> s (cell 1): 1.50649
#> nu (cell 1): 0.428332
#> theta (batch 1): 0.451402
#> betas:    0.7658   0.0333   0.8473   0.5574   0.6939   0.2202   0.2971   0.0721  -0.0913  -0.6681  -0.5752  -0.4598
#> 
#> sigma: 0.126345
#> lambda (gene 1): 0.539132
#> epsilon (gene 1): 0.263867
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 689.134
#> delta (gene 1): 0.667208
#> phi (cell 1): 0.760479
#> s (cell 1): 0.71282
#> nu (cell 1): 0.428368
#> theta (batch 1): 0.602746
#> betas:    0.9149  -0.0539   0.6423   0.5908   0.6365   0.5611  -0.1448   0.5025  -0.1275  -0.3954  -0.3991  -0.1477
#> 
#> sigma: 0.157731
#> lambda (gene 1): 0.696598
#> epsilon (gene 1): -0.239659
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 626.43
#> delta (gene 1): 0.688887
#> phi (cell 1): 0.848309
#> s (cell 1): 0.655156
#> nu (cell 1): 0.433634
#> theta (batch 1): 0.56329
#> betas:    1.0047  -0.1792  -0.0359   1.1307   0.8720   0.2656   0.3678   0.1138   0.4342  -0.1459  -0.1780   0.1049
#> 
#> sigma: 0.155151
#> lambda (gene 1): 1.11468
#> epsilon (gene 1): -0.253858
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 639.945
#> delta (gene 1): 0.700641
#> phi (cell 1): 0.797385
#> s (cell 1): 0.489788
#> nu (cell 1): 0.41732
#> theta (batch 1): 0.197498
#> betas:    0.5029  -0.0117   0.6804   0.0992   1.2353   0.4326  -0.0014   0.4421  -0.2040  -0.0760  -0.5833  -0.3245
#> 
#> sigma: 0.129658
#> lambda (gene 1): 0.501365
#> epsilon (gene 1): -0.127847
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 672.377
#> delta (gene 1): 0.750657
#> phi (cell 1): 0.849095
#> s (cell 1): 0.51507
#> nu (cell 1): 0.434254
#> theta (batch 1): 0.359382
#> betas:    1.0859  -0.1560   0.6243   0.7936   0.8235   0.3020   0.7448  -0.2233   0.0083   0.2159  -0.4004  -0.1870
#> 
#> sigma: 0.179968
#> lambda (gene 1): 2.05087
#> epsilon (gene 1): -0.174229
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.774
#> delta (gene 1): 0.677364
#> phi (cell 1): 0.835354
#> s (cell 1): 0.404145
#> nu (cell 1): 0.409775
#> theta (batch 1): 0.0556431
#> betas:    1.1501  -0.1620   0.6566   0.6365   0.6830   0.3663   0.4571   0.6196  -0.2822  -0.0576  -0.6231   0.4170
#> 
#> sigma: 0.150406
#> lambda (gene 1): 1.18076
#> epsilon (gene 1): -0.0153919
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 673.429
#> delta (gene 1): 0.590463
#> phi (cell 1): 0.882094
#> s (cell 1): 0.470445
#> nu (cell 1): 0.425574
#> theta (batch 1): 0.0599162
#> betas:    0.6927  -0.0492   0.6457   0.7307   0.8390   0.4402   0.2420   0.5118  -0.2325  -0.3717   0.0834  -0.6301
#> 
#> sigma: 0.145429
#> lambda (gene 1): 0.565345
#> epsilon (gene 1): -0.368655
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 813.053
#> delta (gene 1): 0.637231
#> phi (cell 1): 0.861862
#> s (cell 1): 0.569517
#> nu (cell 1): 0.422426
#> theta (batch 1): 0.141047
#> betas:    0.4213   0.0118   0.7768   0.8665   0.7762   0.4755   0.0747   0.8238  -0.8405   0.4126  -0.8849  -0.2672
#> 
#> sigma: 0.137437
#> lambda (gene 1): 0.782866
#> epsilon (gene 1): -0.241155
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 666.362
#> delta (gene 1): 0.655573
#> phi (cell 1): 0.794213
#> s (cell 1): 0.220457
#> nu (cell 1): 0.414157
#> theta (batch 1): 0.299312
#> betas:    0.6787  -0.1324   0.5789   0.7403   0.7957   0.5150   0.4535   0.4240  -0.2950   0.3051  -0.2247   0.0711
#> 
#> sigma: 0.183862
#> lambda (gene 1): 2.21922
#> epsilon (gene 1): -0.288533
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 799.308
#> delta (gene 1): 0.458054
#> phi (cell 1): 0.66737
#> s (cell 1): 0.439944
#> nu (cell 1): 0.42671
#> theta (batch 1): 0.351749
#> betas:    0.5067  -0.0135   0.9475   0.8524   0.7478   0.4807   0.4254  -0.1694   0.5775  -0.3813  -0.8315  -0.2191
#> 
#> sigma: 0.16407
#> lambda (gene 1): 0.261691
#> epsilon (gene 1): -0.47497
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 806.495
#> delta (gene 1): 0.773422
#> phi (cell 1): 0.739627
#> s (cell 1): 0.342535
#> nu (cell 1): 0.427076
#> theta (batch 1): 0.560302
#> betas:    1.1069  -0.1148   0.7852   0.6192   0.2301   0.5639   0.1943   0.3083  -0.2651  -0.0789  -0.0847  -0.7175
#> 
#> sigma: 0.185407
#> lambda (gene 1): 0.293934
#> epsilon (gene 1): -0.148832
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 849.031
#> delta (gene 1): 0.622759
#> phi (cell 1): 0.688107
#> s (cell 1): 0.540382
#> nu (cell 1): 0.418133
#> theta (batch 1): 0.303762
#> betas:    0.7784  -0.0832   1.0339   0.8287   0.9701   0.1470   0.8447   0.2334  -0.2981  -0.0599  -0.1357  -0.5227
#> 
#> sigma: 0.169312
#> lambda (gene 1): 1.03264
#> epsilon (gene 1): -0.268514
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 654.89
#> delta (gene 1): 0.668713
#> phi (cell 1): 0.682013
#> s (cell 1): 0.66369
#> nu (cell 1): 0.426845
#> theta (batch 1): 0.467879
#> betas:    0.9888  -0.2070   0.0977   1.3347   0.0863   0.4766   1.0647  -0.0040   0.2757   0.0817  -0.6874   0.8492
#> 
#> sigma: 0.204953
#> lambda (gene 1): 0.710702
#> epsilon (gene 1): -0.0344537
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 507.179
#> delta (gene 1): 0.547497
#> phi (cell 1): 0.742631
#> s (cell 1): 0.273069
#> nu (cell 1): 0.427028
#> theta (batch 1): 0.19999
#> betas:    0.8022  -0.0616   0.7768   1.0405   0.6152   0.1871   0.5715   0.0478  -0.0358  -0.0676  -0.6208  -0.0971
#> 
#> sigma: 0.131882
#> lambda (gene 1): 1.10392
#> epsilon (gene 1): -0.432972
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 583.888
#> delta (gene 1): 0.721321
#> phi (cell 1): 0.719655
#> s (cell 1): 0.552557
#> nu (cell 1): 0.429227
#> theta (batch 1): 0.214991
#> betas:    1.1134  -0.1619   0.4887   0.5637   0.7264   0.3622   0.5039  -0.1741   0.5377  -0.3885  -0.2119  -0.0616
#> 
#> sigma: 0.112102
#> lambda (gene 1): 1.02238
#> epsilon (gene 1): -0.126671
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 831.781
#> delta (gene 1): 0.893169
#> phi (cell 1): 0.734757
#> s (cell 1): 0.367533
#> nu (cell 1): 0.434434
#> theta (batch 1): 0.030457
#> betas:    1.1445  -0.1024   0.4668   0.7057   0.3687   0.8298  -0.0732   0.0195  -0.0067  -0.2816  -0.7752   0.0453
#> 
#> sigma: 0.185424
#> lambda (gene 1): 0.660789
#> epsilon (gene 1): 0.282213
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 620.454
#> delta (gene 1): 0.898542
#> phi (cell 1): 0.76089
#> s (cell 1): 0.407898
#> nu (cell 1): 0.420107
#> theta (batch 1): 0.0167822
#> betas:    0.6359  -0.0772   1.1693   0.4035   0.6625   0.7422   0.0419   0.5094   0.1251  -0.3888  -0.0463  -0.2886
#> 
#> sigma: 0.170489
#> lambda (gene 1): 3.19468
#> epsilon (gene 1): 0.0852944
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 673.148
#> delta (gene 1): 0.562152
#> phi (cell 1): 0.71801
#> s (cell 1): 0.475305
#> nu (cell 1): 0.432935
#> theta (batch 1): 0.028406
#> betas:    0.6414  -0.0119   0.5278   0.8328   0.3355   0.5629   0.1069   0.4227  -0.1400  -0.4861  -0.2441  -0.3412
#> 
#> sigma: 0.172273
#> lambda (gene 1): 0.998251
#> epsilon (gene 1): -0.372387
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 763.456
#> delta (gene 1): 0.625213
#> phi (cell 1): 0.697931
#> s (cell 1): 0.455703
#> nu (cell 1): 0.425052
#> theta (batch 1): 0.0227239
#> betas:    1.2009  -0.1219   0.5658   0.7356   0.7205   0.2351   0.3967  -0.1105   0.1693  -0.1600  -0.5403  -0.2199
#> 
#> sigma: 0.149211
#> lambda (gene 1): 0.692305
#> epsilon (gene 1): -0.324562
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 684.325
#> delta (gene 1): 1.07732
#> phi (cell 1): 0.780727
#> s (cell 1): 0.352166
#> nu (cell 1): 0.429111
#> theta (batch 1): 0.062925
#> betas:    0.8636  -0.0642   0.6160   0.7445   0.5623   0.5675   0.1905   0.0259   0.3494  -0.7042  -0.3136   0.1831
#> 
#> sigma: 0.126506
#> lambda (gene 1): 0.521107
#> epsilon (gene 1): 0.284457
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 563.872
#> delta (gene 1): 0.704504
#> phi (cell 1): 0.819995
#> s (cell 1): 0.640556
#> nu (cell 1): 0.43266
#> theta (batch 1): 0.153827
#> betas:    1.2954  -0.1399   0.4794   0.8429   0.2976   0.6486   0.0569   0.2665  -0.5090   0.1831  -0.4771  -0.3371
#> 
#> sigma: 0.201901
#> lambda (gene 1): 0.808863
#> epsilon (gene 1): -0.164907
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 9980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 780.735
#> delta (gene 1): 0.741852
#> phi (cell 1): 1.00692
#> s (cell 1): 0.380964
#> nu (cell 1): 0.407835
#> theta (batch 1): 0.091183
#> betas:    1.3087  -0.1448   0.6820   0.8364   0.3681   0.5213   0.3900   0.0344  -0.2106  -0.0628  -0.0565  -0.5540
#> 
#> sigma: 0.146824
#> lambda (gene 1): 1.26242
#> epsilon (gene 1): -0.262341
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 790.452
#> delta (gene 1): 0.764193
#> phi (cell 1): 0.914675
#> s (cell 1): 0.417613
#> nu (cell 1): 0.417065
#> theta (batch 1): 0.216225
#> betas:    0.5848  -0.0827   0.4016   0.9020   0.9122   0.4752   0.5533   0.3867  -0.4653   0.2441  -0.1724  -0.4685
#> 
#> sigma: 0.111204
#> lambda (gene 1): 0.793792
#> epsilon (gene 1): -0.0782904
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 654.314
#> delta (gene 1): 0.816177
#> phi (cell 1): 0.844309
#> s (cell 1): 0.248894
#> nu (cell 1): 0.4097
#> theta (batch 1): 0.471181
#> betas:    0.8759  -0.0409   1.0533   0.4320   0.1772   0.9969  -0.1016   0.0082  -0.0157  -0.2930  -0.2990  -0.6107
#> 
#> sigma: 0.124621
#> lambda (gene 1): 0.696721
#> epsilon (gene 1): -0.0999011
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 682.915
#> delta (gene 1): 0.535171
#> phi (cell 1): 0.780752
#> s (cell 1): 0.455255
#> nu (cell 1): 0.433679
#> theta (batch 1): 0.459058
#> betas:    0.8144   0.0528   1.2919   0.4733   0.6612   0.1154   0.2196   0.0208  -0.4008  -0.0107  -1.1352  -0.6415
#> 
#> sigma: 0.181565
#> lambda (gene 1): 1.56266
#> epsilon (gene 1): -0.479611
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 598.384
#> delta (gene 1): 0.475894
#> phi (cell 1): 0.808185
#> s (cell 1): 0.593641
#> nu (cell 1): 0.415616
#> theta (batch 1): 0.135948
#> betas:    0.8731   0.0140   1.1336   0.0786   0.6603   0.7358   0.2926  -0.5614   0.2422  -0.5995  -0.4922  -0.6564
#> 
#> sigma: 0.201754
#> lambda (gene 1): 2.28423
#> epsilon (gene 1): -0.596075
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 812.221
#> delta (gene 1): 0.536016
#> phi (cell 1): 0.782779
#> s (cell 1): 0.428256
#> nu (cell 1): 0.423679
#> theta (batch 1): 0.0239365
#> betas:    1.2582  -0.2117   0.5636   0.8470   0.2824   0.3922   0.8270  -0.0453   0.2362  -0.0669  -0.0520  -0.2230
#> 
#> sigma: 0.131912
#> lambda (gene 1): 0.522731
#> epsilon (gene 1): -0.474715
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 740.985
#> delta (gene 1): 0.722512
#> phi (cell 1): 0.750293
#> s (cell 1): 0.512523
#> nu (cell 1): 0.441658
#> theta (batch 1): 0.0473849
#> betas:    0.7021  -0.0516   0.1857   1.0562   0.4224   0.6136   0.3709  -0.0330   0.0819  -0.1811  -0.5616   0.1063
#> 
#> sigma: 0.157448
#> lambda (gene 1): 0.593728
#> epsilon (gene 1): -0.167071
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 775.867
#> delta (gene 1): 1.12611
#> phi (cell 1): 0.742347
#> s (cell 1): 0.481949
#> nu (cell 1): 0.406004
#> theta (batch 1): 0.02209
#> betas:    0.8618  -0.1399   0.8683   0.6530   0.6931   0.6258   0.1313   0.3237   0.3264  -0.2691  -0.2570  -0.0208
#> 
#> sigma: 0.171816
#> lambda (gene 1): 2.02474
#> epsilon (gene 1): 0.371318
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 743.299
#> delta (gene 1): 0.961488
#> phi (cell 1): 0.773991
#> s (cell 1): 0.480108
#> nu (cell 1): 0.422626
#> theta (batch 1): 0.0195547
#> betas:    1.0437  -0.1193   0.7140   0.7897   0.7529   0.1356   0.3975   0.0427   0.3199  -0.6603  -0.0809  -0.2396
#> 
#> sigma: 0.166076
#> lambda (gene 1): 1.53045
#> epsilon (gene 1): 0.269414
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 783.586
#> delta (gene 1): 1.20965
#> phi (cell 1): 0.784003
#> s (cell 1): 0.500543
#> nu (cell 1): 0.417663
#> theta (batch 1): 0.0192725
#> betas:    0.6514   0.0047   0.5078   0.7888   0.8299   0.3527   0.2345   0.2429  -0.2916   0.0536  -0.8302  -0.4992
#> 
#> sigma: 0.164168
#> lambda (gene 1): 0.808076
#> epsilon (gene 1): 0.36166
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 767.58
#> delta (gene 1): 0.734118
#> phi (cell 1): 0.806875
#> s (cell 1): 0.420464
#> nu (cell 1): 0.421205
#> theta (batch 1): 0.00696614
#> betas:    0.9566  -0.0481   0.5956   0.5154   1.1024   0.5493  -0.5316   0.3236   0.1330  -0.1698  -0.8562  -0.1180
#> 
#> sigma: 0.165901
#> lambda (gene 1): 1.83556
#> epsilon (gene 1): -0.227389
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 654.601
#> delta (gene 1): 0.690206
#> phi (cell 1): 0.772282
#> s (cell 1): 0.430559
#> nu (cell 1): 0.428073
#> theta (batch 1): 0.00595331
#> betas:    1.0052  -0.1139   0.0128   0.9136   0.6172   0.4543   0.3004   0.1427   0.0260  -0.5805  -0.2176   0.0901
#> 
#> sigma: 0.196048
#> lambda (gene 1): 2.37874
#> epsilon (gene 1): 0.0408195
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 874.115
#> delta (gene 1): 0.610317
#> phi (cell 1): 0.938971
#> s (cell 1): 0.388504
#> nu (cell 1): 0.412092
#> theta (batch 1): 0.0198386
#> betas:    1.6379  -0.1617   0.6112   0.3553   0.4823   0.3600  -0.0328   0.2628  -0.0485  -0.6572   0.1089  -0.6277
#> 
#> sigma: 0.152755
#> lambda (gene 1): 1.70856
#> epsilon (gene 1): -0.32566
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 779.161
#> delta (gene 1): 0.849924
#> phi (cell 1): 0.934583
#> s (cell 1): 0.523674
#> nu (cell 1): 0.413629
#> theta (batch 1): 0.035453
#> betas:    0.6577  -0.0940   0.3467   1.2351   0.6884   0.6461   0.2568   0.3484  -0.0201  -0.2568   0.1523  -0.2952
#> 
#> sigma: 0.149928
#> lambda (gene 1): 1.00466
#> epsilon (gene 1): -0.040221
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 599.787
#> delta (gene 1): 0.825616
#> phi (cell 1): 0.871452
#> s (cell 1): 0.818011
#> nu (cell 1): 0.42209
#> theta (batch 1): 0.10268
#> betas:    0.5826  -0.0811   0.8170   0.8538   0.1669   1.4514   0.2177   0.2764  -0.0041  -0.3710   0.2185  -0.4740
#> 
#> sigma: 0.157123
#> lambda (gene 1): 0.724407
#> epsilon (gene 1): 0.0389586
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 642.674
#> delta (gene 1): 0.842299
#> phi (cell 1): 0.819328
#> s (cell 1): 0.127805
#> nu (cell 1): 0.401084
#> theta (batch 1): 0.359317
#> betas:    0.9079  -0.1000   0.2892   1.0695   0.5217   0.2667   0.8685   0.1314   0.0496  -0.4171  -0.6316   0.5326
#> 
#> sigma: 0.246177
#> lambda (gene 1): 0.795315
#> epsilon (gene 1): 0.236797
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 762.81
#> delta (gene 1): 0.82078
#> phi (cell 1): 0.791783
#> s (cell 1): 0.650663
#> nu (cell 1): 0.431436
#> theta (batch 1): 0.289413
#> betas:    0.4149  -0.0112   0.6230   0.9083   1.1430   0.3056   0.4430   0.1439   0.1630  -0.3178  -0.5613   0.2308
#> 
#> sigma: 0.154132
#> lambda (gene 1): 1.05878
#> epsilon (gene 1): -0.0130477
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 813.024
#> delta (gene 1): 0.525906
#> phi (cell 1): 0.855719
#> s (cell 1): 0.407573
#> nu (cell 1): 0.421944
#> theta (batch 1): 0.387696
#> betas:    0.6344  -0.0771   0.5433   0.7369   1.0502   0.3613   0.5505   0.0403   0.2204  -0.3299  -0.2984   0.0981
#> 
#> sigma: 0.193933
#> lambda (gene 1): 1.04037
#> epsilon (gene 1): -0.408443
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 648.924
#> delta (gene 1): 0.861038
#> phi (cell 1): 0.696497
#> s (cell 1): 0.329705
#> nu (cell 1): 0.423348
#> theta (batch 1): 0.505761
#> betas:    0.4449  -0.0881   0.5005   1.1672   0.8072   0.4431   0.7527   0.0452   0.2325  -0.3273  -0.1516   0.0694
#> 
#> sigma: 0.140908
#> lambda (gene 1): 0.722549
#> epsilon (gene 1): 0.225256
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 573.049
#> delta (gene 1): 0.490324
#> phi (cell 1): 0.713523
#> s (cell 1): 1.61908
#> nu (cell 1): 0.427428
#> theta (batch 1): 0.634712
#> betas:    0.9161  -0.1349   0.6222   0.8922   0.7430   0.3071   0.3313   0.2645  -0.2280   0.0582  -0.2372   0.0924
#> 
#> sigma: 0.123253
#> lambda (gene 1): 0.520831
#> epsilon (gene 1): -0.555289
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 722.19
#> delta (gene 1): 0.568415
#> phi (cell 1): 0.732243
#> s (cell 1): 0.2593
#> nu (cell 1): 0.432642
#> theta (batch 1): 0.160092
#> betas:    1.2392  -0.2008   0.3004   0.6350   0.5251   0.5893   0.8974  -0.3063   0.2189  -0.0101   0.2115  -0.6006
#> 
#> sigma: 0.19299
#> lambda (gene 1): 0.892079
#> epsilon (gene 1): -0.5536
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 877.046
#> delta (gene 1): 0.633715
#> phi (cell 1): 0.750319
#> s (cell 1): 0.717102
#> nu (cell 1): 0.439559
#> theta (batch 1): 0.166032
#> betas:    0.9582  -0.1669   0.8200   0.5709   0.9490   0.0895   0.7662   0.6452  -0.5063   0.1370   0.2533  -0.2999
#> 
#> sigma: 0.148278
#> lambda (gene 1): 0.980554
#> epsilon (gene 1): -0.380169
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 647.727
#> delta (gene 1): 0.575807
#> phi (cell 1): 0.685566
#> s (cell 1): 0.276035
#> nu (cell 1): 0.429768
#> theta (batch 1): 0.185726
#> betas:    0.6663  -0.0109  -0.1859   0.7656   0.8309   0.2864   0.5149  -0.2792   0.0691  -0.5922  -0.1819  -0.5319
#> 
#> sigma: 0.146914
#> lambda (gene 1): 3.22976
#> epsilon (gene 1): -0.267143
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 811.633
#> delta (gene 1): 0.656811
#> phi (cell 1): 0.707799
#> s (cell 1): 0.481622
#> nu (cell 1): 0.43858
#> theta (batch 1): 0.17729
#> betas:    1.0785  -0.1187   0.0475   1.0812   0.8677   0.2888   0.2159   0.1480   0.2747  -0.4838  -0.4556   0.2283
#> 
#> sigma: 0.124525
#> lambda (gene 1): 0.807635
#> epsilon (gene 1): -0.164576
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 598.518
#> delta (gene 1): 0.788885
#> phi (cell 1): 0.636379
#> s (cell 1): 0.315741
#> nu (cell 1): 0.424818
#> theta (batch 1): 0.25958
#> betas:    0.9800  -0.0063   1.4070   0.3380   0.9688   0.0125   0.5300  -0.2773   0.2561  -0.6626  -0.6861  -0.4345
#> 
#> sigma: 0.186515
#> lambda (gene 1): 1.66455
#> epsilon (gene 1): 0.00819745
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 703.613
#> delta (gene 1): 0.791803
#> phi (cell 1): 0.582917
#> s (cell 1): 0.519826
#> nu (cell 1): 0.417895
#> theta (batch 1): 0.104987
#> betas:    0.4849   0.0700   1.0043   0.8012   0.3282   0.4914   0.4302  -0.2008  -0.1535  -0.3275  -0.9000  -0.6189
#> 
#> sigma: 0.134828
#> lambda (gene 1): 1.72287
#> epsilon (gene 1): 0.134633
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 815.771
#> delta (gene 1): 0.987252
#> phi (cell 1): 0.626645
#> s (cell 1): 0.926328
#> nu (cell 1): 0.422306
#> theta (batch 1): 0.178253
#> betas:    0.5766  -0.1144   0.7202   0.6414   0.8174   0.7318   0.2101   0.3105   0.3286  -0.4617   0.2090  -0.1534
#> 
#> sigma: 0.179902
#> lambda (gene 1): 1.33502
#> epsilon (gene 1): 0.228466
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 746.646
#> delta (gene 1): 0.636731
#> phi (cell 1): 0.631146
#> s (cell 1): 0.408906
#> nu (cell 1): 0.422959
#> theta (batch 1): 0.370573
#> betas:    0.5317   0.0049   0.8689   0.6392   1.1597   0.4370  -0.0412   0.5357  -0.1280  -0.4264  -0.6556  -0.3770
#> 
#> sigma: 0.17201
#> lambda (gene 1): 0.227687
#> epsilon (gene 1): -0.0114955
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.062
#> delta (gene 1): 0.983371
#> phi (cell 1): 0.662058
#> s (cell 1): 0.667473
#> nu (cell 1): 0.41961
#> theta (batch 1): 0.423574
#> betas:    1.0748  -0.2260   0.3324   1.0201   1.0160   0.5672   0.4847   0.1267   0.2070   0.3290  -0.4833   0.2934
#> 
#> sigma: 0.142173
#> lambda (gene 1): 0.622716
#> epsilon (gene 1): 0.188113
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 607.751
#> delta (gene 1): 0.725706
#> phi (cell 1): 0.68076
#> s (cell 1): 0.354002
#> nu (cell 1): 0.411307
#> theta (batch 1): 0.518119
#> betas:    0.9865  -0.0882  -0.0982   0.6051   0.9791   0.5620   0.0732  -0.0122   0.1445  -0.4582  -0.5084  -0.1329
#> 
#> sigma: 0.152889
#> lambda (gene 1): 1.05949
#> epsilon (gene 1): 0.0672078
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 895.723
#> delta (gene 1): 0.625766
#> phi (cell 1): 0.763142
#> s (cell 1): 0.488965
#> nu (cell 1): 0.420763
#> theta (batch 1): 0.622666
#> betas:    0.5288  -0.0405   0.0844   1.1473   0.7487   0.8259   0.1683   0.0240   0.1071  -0.1865  -0.1692  -0.6638
#> 
#> sigma: 0.159945
#> lambda (gene 1): 0.730223
#> epsilon (gene 1): -0.30262
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 603.103
#> delta (gene 1): 0.526965
#> phi (cell 1): 0.745206
#> s (cell 1): 0.793521
#> nu (cell 1): 0.434301
#> theta (batch 1): 0.350357
#> betas:    1.3097  -0.0814   0.5287   0.5157   0.4546   0.4446  -0.2279   0.3100  -0.1085  -0.5926  -0.3863  -0.5433
#> 
#> sigma: 0.140139
#> lambda (gene 1): 0.880709
#> epsilon (gene 1): -0.354509
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 758.046
#> delta (gene 1): 0.608296
#> phi (cell 1): 0.752168
#> s (cell 1): 0.628211
#> nu (cell 1): 0.42817
#> theta (batch 1): 0.330483
#> betas:    0.6911  -0.1900   0.7766   1.0808   0.5740   0.8490   0.6457   0.4074  -0.0008  -0.0797   0.2039   0.3054
#> 
#> sigma: 0.177732
#> lambda (gene 1): 1.18622
#> epsilon (gene 1): -0.199334
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 550.867
#> delta (gene 1): 0.592218
#> phi (cell 1): 0.740697
#> s (cell 1): 0.461015
#> nu (cell 1): 0.426738
#> theta (batch 1): 0.157299
#> betas:    0.6988  -0.0343   0.4929   0.9048   1.0333   0.0363   0.6735  -0.0965   0.2750  -0.4523  -0.2916  -0.5855
#> 
#> sigma: 0.1578
#> lambda (gene 1): 0.45543
#> epsilon (gene 1): -0.332177
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 712.299
#> delta (gene 1): 0.497912
#> phi (cell 1): 0.616434
#> s (cell 1): 0.476399
#> nu (cell 1): 0.431701
#> theta (batch 1): 0.230578
#> betas:    0.9302  -0.0238   0.5442   0.7843   0.8090   0.0177   0.4430  -0.0468  -0.2804  -0.1841  -0.3223  -0.8116
#> 
#> sigma: 0.120899
#> lambda (gene 1): 0.625435
#> epsilon (gene 1): -0.627273
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 709.615
#> delta (gene 1): 0.672814
#> phi (cell 1): 0.632656
#> s (cell 1): 0.35415
#> nu (cell 1): 0.416077
#> theta (batch 1): 0.127587
#> betas:    0.4016  -0.0470   1.0332   0.5780   0.8971   1.0731  -0.4569   0.7406   0.1397  -0.4347  -0.2619  -0.0509
#> 
#> sigma: 0.200125
#> lambda (gene 1): 1.07794
#> epsilon (gene 1): -0.0725372
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 589.115
#> delta (gene 1): 0.4707
#> phi (cell 1): 0.636801
#> s (cell 1): 0.45681
#> nu (cell 1): 0.433302
#> theta (batch 1): 0.074994
#> betas:    0.9626  -0.0443   1.0581   0.5478   0.3481   0.2601   0.7141  -0.0206  -0.7670   0.2395  -0.6353  -0.3987
#> 
#> sigma: 0.166991
#> lambda (gene 1): 0.296881
#> epsilon (gene 1): -0.582291
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 767.205
#> delta (gene 1): 0.669876
#> phi (cell 1): 0.64223
#> s (cell 1): 0.890197
#> nu (cell 1): 0.421376
#> theta (batch 1): 0.196308
#> betas:    0.9394  -0.1192   0.4066   0.9774   0.7086   0.5959  -0.0103   0.3746  -0.1348   0.1789  -0.2646  -0.4744
#> 
#> sigma: 0.226329
#> lambda (gene 1): 0.839527
#> epsilon (gene 1): -0.377549
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 639.501
#> delta (gene 1): 0.695264
#> phi (cell 1): 0.673476
#> s (cell 1): 1.36824
#> nu (cell 1): 0.430909
#> theta (batch 1): 0.39644
#> betas:    1.0356  -0.1028   0.4668   0.6839   0.5417   0.4327   0.6424  -0.3519   0.3247  -0.7046   0.0859  -0.5326
#> 
#> sigma: 0.220448
#> lambda (gene 1): 0.732264
#> epsilon (gene 1): -0.0837071
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 710.546
#> delta (gene 1): 0.690023
#> phi (cell 1): 0.731124
#> s (cell 1): 1.48537
#> nu (cell 1): 0.432252
#> theta (batch 1): 0.165124
#> betas:    0.8286  -0.1072   0.3980   0.8342   0.8480   0.2589   0.4041   0.3319  -0.2922   0.4580  -0.6953  -0.0491
#> 
#> sigma: 0.164188
#> lambda (gene 1): 0.325482
#> epsilon (gene 1): -0.293063
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 630.102
#> delta (gene 1): 0.659101
#> phi (cell 1): 0.540395
#> s (cell 1): 0.640474
#> nu (cell 1): 0.413759
#> theta (batch 1): 0.290302
#> betas:    1.0035  -0.1342   0.6679   0.4373   0.8282   0.3666   0.6698  -0.1722   0.3288  -0.4012  -0.1601  -0.0856
#> 
#> sigma: 0.150883
#> lambda (gene 1): 0.904323
#> epsilon (gene 1): -0.182242
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.829
#> delta (gene 1): 0.706526
#> phi (cell 1): 0.547993
#> s (cell 1): 0.441051
#> nu (cell 1): 0.425233
#> theta (batch 1): 0.074919
#> betas:    0.6324  -0.0720   0.0925   1.3242   0.6539   0.6260   0.2551  -0.0847   0.1724  -0.0250  -0.6386   0.0851
#> 
#> sigma: 0.146842
#> lambda (gene 1): 2.23651
#> epsilon (gene 1): -0.117978
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 612.716
#> delta (gene 1): 0.48078
#> phi (cell 1): 0.613056
#> s (cell 1): 0.677964
#> nu (cell 1): 0.409379
#> theta (batch 1): 0.133786
#> betas:    1.1679  -0.0495   0.3915   0.3521   0.6051   0.2776   0.0717   0.4911  -0.6075  -0.5074  -0.1066  -0.9544
#> 
#> sigma: 0.155651
#> lambda (gene 1): 1.5859
#> epsilon (gene 1): -0.418995
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 543.628
#> delta (gene 1): 0.909219
#> phi (cell 1): 0.713724
#> s (cell 1): 0.781137
#> nu (cell 1): 0.426118
#> theta (batch 1): 0.215684
#> betas:    0.7785  -0.0643   0.7847   0.6919   0.6999   0.3528   0.4127   0.1266  -0.3755   0.1692  -0.5937  -0.2120
#> 
#> sigma: 0.159437
#> lambda (gene 1): 0.823669
#> epsilon (gene 1): 0.101743
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 605.938
#> delta (gene 1): 1.33519
#> phi (cell 1): 0.698228
#> s (cell 1): 0.368648
#> nu (cell 1): 0.434273
#> theta (batch 1): 0.0387149
#> betas:    0.7865  -0.0176   0.4843   0.7964   0.3869   0.7990  -0.0688   0.3054  -0.5870   0.2099  -0.7212  -0.5835
#> 
#> sigma: 0.137406
#> lambda (gene 1): 0.901998
#> epsilon (gene 1): 0.447123
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 693.475
#> delta (gene 1): 0.546524
#> phi (cell 1): 0.731285
#> s (cell 1): 1.26522
#> nu (cell 1): 0.424898
#> theta (batch 1): 0.120353
#> betas:    0.5212  -0.0560   0.9434   1.0502   0.5002   0.6482   0.5751  -0.3045   0.6282  -0.4560  -0.4266   0.1244
#> 
#> sigma: 0.111814
#> lambda (gene 1): 1.73749
#> epsilon (gene 1): -0.374332
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 654.535
#> delta (gene 1): 0.66824
#> phi (cell 1): 0.682191
#> s (cell 1): 0.350378
#> nu (cell 1): 0.418136
#> theta (batch 1): 0.0436028
#> betas:    0.7583  -0.0084   0.5640   0.2741   1.0107   0.4037   0.1986   0.0396   0.1992  -0.6741  -0.0793  -0.5303
#> 
#> sigma: 0.163876
#> lambda (gene 1): 1.04797
#> epsilon (gene 1): -0.365161
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 647.227
#> delta (gene 1): 0.581217
#> phi (cell 1): 0.675853
#> s (cell 1): 0.484948
#> nu (cell 1): 0.415091
#> theta (batch 1): 0.0211111
#> betas:    1.2162  -0.1671   0.7602   0.5084   0.7953  -0.2758   0.6387   0.2342  -0.4147   0.2364  -0.4592   0.2581
#> 
#> sigma: 0.142854
#> lambda (gene 1): 0.724029
#> epsilon (gene 1): -0.409256
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 822.822
#> delta (gene 1): 0.737138
#> phi (cell 1): 0.600741
#> s (cell 1): 0.409275
#> nu (cell 1): 0.431397
#> theta (batch 1): 0.00348931
#> betas:    1.0007   0.0400   0.4862   0.4101   0.6113   0.3604   0.0886  -0.1385  -0.3766  -0.5485  -0.7903  -0.5655
#> 
#> sigma: 0.176369
#> lambda (gene 1): 0.745895
#> epsilon (gene 1): 0.0115438
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 767.512
#> delta (gene 1): 0.537773
#> phi (cell 1): 0.670922
#> s (cell 1): 0.513756
#> nu (cell 1): 0.412993
#> theta (batch 1): 0.021896
#> betas:    1.3074  -0.2185   0.3993   1.1849   0.1447   0.6069   0.6406   0.0141   0.0016   0.0930  -0.5180   0.3387
#> 
#> sigma: 0.138187
#> lambda (gene 1): 1.52773
#> epsilon (gene 1): -0.29761
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 10980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 751.235
#> delta (gene 1): 0.902502
#> phi (cell 1): 0.699374
#> s (cell 1): 0.312738
#> nu (cell 1): 0.426432
#> theta (batch 1): 0.0767168
#> betas:    0.7781  -0.1329   0.4315   0.9746   0.7992   0.6293  -0.0041   0.5043   0.2666  -0.1691  -0.4053   0.0467
#> 
#> sigma: 0.154355
#> lambda (gene 1): 1.91257
#> epsilon (gene 1): 0.186527
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 812.146
#> delta (gene 1): 0.519361
#> phi (cell 1): 0.682134
#> s (cell 1): 0.535508
#> nu (cell 1): 0.412947
#> theta (batch 1): 0.10337
#> betas:    0.7408  -0.0009   0.8708   0.7747   0.7855   0.4800   0.2729   0.1030  -0.3129  -0.2877  -0.6147  -0.5257
#> 
#> sigma: 0.124658
#> lambda (gene 1): 0.126729
#> epsilon (gene 1): -0.311084
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 606.488
#> delta (gene 1): 0.584115
#> phi (cell 1): 0.649902
#> s (cell 1): 0.463059
#> nu (cell 1): 0.432211
#> theta (batch 1): 0.0163583
#> betas:    0.3657  -0.0071   0.9291   0.6845   0.3348   1.1811  -0.1600   0.3835   0.0269  -0.4036  -0.0244  -0.4674
#> 
#> sigma: 0.168687
#> lambda (gene 1): 1.33293
#> epsilon (gene 1): -0.395458
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 895.583
#> delta (gene 1): 0.701365
#> phi (cell 1): 0.687602
#> s (cell 1): 0.402428
#> nu (cell 1): 0.434356
#> theta (batch 1): 0.0604399
#> betas:    1.2991  -0.1579   0.5226   0.7460   0.5984   0.8461   0.0897  -0.0653   0.1471  -0.3101   0.0775  -0.4960
#> 
#> sigma: 0.189002
#> lambda (gene 1): 0.696268
#> epsilon (gene 1): -0.295677
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 791.744
#> delta (gene 1): 0.699853
#> phi (cell 1): 0.703045
#> s (cell 1): 0.518305
#> nu (cell 1): 0.421485
#> theta (batch 1): 0.0719556
#> betas:    1.0192  -0.1496   0.7388   0.7260   0.3827   0.8269   0.2368   0.2170   0.1011  -0.0374  -0.6469   0.3628
#> 
#> sigma: 0.153904
#> lambda (gene 1): 0.5866
#> epsilon (gene 1): -0.0780081
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 733.501
#> delta (gene 1): 0.846938
#> phi (cell 1): 0.843012
#> s (cell 1): 0.44368
#> nu (cell 1): 0.424456
#> theta (batch 1): 0.244388
#> betas:    1.3497  -0.1624  -0.8072   0.9164   0.9465   0.1205   0.2235   0.3657  -0.2167  -0.1871  -0.0193  -0.5520
#> 
#> sigma: 0.216036
#> lambda (gene 1): 2.15258
#> epsilon (gene 1): -0.0178038
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 743.36
#> delta (gene 1): 0.964953
#> phi (cell 1): 0.735006
#> s (cell 1): 0.57156
#> nu (cell 1): 0.417938
#> theta (batch 1): 0.194573
#> betas:    0.9326  -0.1382   0.8074   0.6446   0.6352   0.7213   0.5315  -0.0904   0.1360  -0.1212  -0.2448   0.1348
#> 
#> sigma: 0.118135
#> lambda (gene 1): 1.18983
#> epsilon (gene 1): 0.129877
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 721.936
#> delta (gene 1): 0.708673
#> phi (cell 1): 0.73813
#> s (cell 1): 0.548949
#> nu (cell 1): 0.419975
#> theta (batch 1): 0.378802
#> betas:    1.4862  -0.1425   0.7189   1.3187  -0.0277   0.2675   0.6417  -0.2712  -0.1637  -0.1946  -0.1703  -0.8855
#> 
#> sigma: 0.197811
#> lambda (gene 1): 2.08899
#> epsilon (gene 1): -0.159633
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 597.284
#> delta (gene 1): 0.971562
#> phi (cell 1): 0.779863
#> s (cell 1): 0.346244
#> nu (cell 1): 0.425835
#> theta (batch 1): 0.206203
#> betas:    0.8897  -0.1449   0.8697   0.5604   0.1015   0.8865   0.4417   0.1086   0.1231  -0.1986   0.0345  -0.1034
#> 
#> sigma: 0.138265
#> lambda (gene 1): 1.76149
#> epsilon (gene 1): 0.126723
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 832.105
#> delta (gene 1): 0.615407
#> phi (cell 1): 0.773143
#> s (cell 1): 0.465731
#> nu (cell 1): 0.420359
#> theta (batch 1): 0.123049
#> betas:    1.0680  -0.1034   0.3584   0.7696   0.4502   0.6500   0.0918  -0.1109   0.5353  -0.4715  -0.6833   0.3162
#> 
#> sigma: 0.219547
#> lambda (gene 1): 0.168158
#> epsilon (gene 1): -0.296308
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 568.541
#> delta (gene 1): 0.594663
#> phi (cell 1): 0.787318
#> s (cell 1): 0.55386
#> nu (cell 1): 0.405618
#> theta (batch 1): 0.0488363
#> betas:    0.4637  -0.0917   0.7150   1.1219   1.1487   0.1875   0.7129   0.0031   0.1276   0.0236   0.2339  -0.6085
#> 
#> sigma: 0.186003
#> lambda (gene 1): 0.2475
#> epsilon (gene 1): -0.505873
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 661.92
#> delta (gene 1): 0.591767
#> phi (cell 1): 0.734238
#> s (cell 1): 0.346933
#> nu (cell 1): 0.432076
#> theta (batch 1): 0.0147485
#> betas:    1.0191  -0.1229   0.1172   1.1006   0.7040   0.7189  -0.1058   0.2664  -0.0111   0.1224  -0.6703   0.2119
#> 
#> sigma: 0.172692
#> lambda (gene 1): 0.749614
#> epsilon (gene 1): -0.445934
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 821.125
#> delta (gene 1): 0.619685
#> phi (cell 1): 0.691581
#> s (cell 1): 0.575665
#> nu (cell 1): 0.43067
#> theta (batch 1): 0.131356
#> betas:    0.5004  -0.0094   0.8888   1.0431   0.0496   0.7289   0.5860   0.0606  -0.5872   0.5289  -1.0842  -0.0763
#> 
#> sigma: 0.156703
#> lambda (gene 1): 0.978823
#> epsilon (gene 1): -0.247756
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 642.93
#> delta (gene 1): 0.605869
#> phi (cell 1): 0.716313
#> s (cell 1): 0.529632
#> nu (cell 1): 0.425318
#> theta (batch 1): 0.327211
#> betas:    0.7318  -0.1968   1.0505   0.4014   0.9186   0.5587   0.3273   0.4002   0.3711   0.4549  -0.5456   0.4850
#> 
#> sigma: 0.158438
#> lambda (gene 1): 0.35691
#> epsilon (gene 1): -0.439557
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 753.014
#> delta (gene 1): 0.651256
#> phi (cell 1): 0.830301
#> s (cell 1): 0.744447
#> nu (cell 1): 0.415239
#> theta (batch 1): 0.241309
#> betas:    0.9704  -0.1925   0.6566   0.9843   0.6472   0.3823   0.7428  -0.0042   0.1032   0.0605  -0.2718   0.3574
#> 
#> sigma: 0.18155
#> lambda (gene 1): 1.24869
#> epsilon (gene 1): -0.167908
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 551.112
#> delta (gene 1): 0.547434
#> phi (cell 1): 0.896931
#> s (cell 1): 1.09807
#> nu (cell 1): 0.420361
#> theta (batch 1): 0.373086
#> betas:    0.8249  -0.1382   0.5189   0.7814   0.5278   0.6024   0.6142  -0.1485  -0.0910   0.2571  -0.4589   0.1921
#> 
#> sigma: 0.19922
#> lambda (gene 1): 1.18328
#> epsilon (gene 1): -0.402547
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 677.941
#> delta (gene 1): 0.925187
#> phi (cell 1): 0.82377
#> s (cell 1): 0.76659
#> nu (cell 1): 0.425607
#> theta (batch 1): 0.369756
#> betas:    1.0290  -0.1189   0.4851   0.9417   0.3321   0.8603  -0.0428   0.1958   0.4403  -0.2604  -0.7575   0.2317
#> 
#> sigma: 0.165686
#> lambda (gene 1): 1.33794
#> epsilon (gene 1): 0.125583
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 584.523
#> delta (gene 1): 0.847212
#> phi (cell 1): 0.779523
#> s (cell 1): 0.499277
#> nu (cell 1): 0.43522
#> theta (batch 1): 0.260518
#> betas:    1.1786  -0.1398   0.4186   0.3820   0.8515   0.0214   0.5671   0.4374  -0.5835   0.0586  -0.2076  -0.1784
#> 
#> sigma: 0.190182
#> lambda (gene 1): 0.45603
#> epsilon (gene 1): 0.00413914
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 615.857
#> delta (gene 1): 0.71251
#> phi (cell 1): 0.709899
#> s (cell 1): 0.405576
#> nu (cell 1): 0.430124
#> theta (batch 1): 0.226081
#> betas:    1.3777  -0.0015   0.7757   0.5068   0.6219   0.0086   0.1951  -0.3717  -0.1817  -0.2413  -1.1430  -0.6159
#> 
#> sigma: 0.135823
#> lambda (gene 1): 0.742123
#> epsilon (gene 1): -0.22537
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 701.545
#> delta (gene 1): 0.664423
#> phi (cell 1): 0.734138
#> s (cell 1): 0.511083
#> nu (cell 1): 0.425004
#> theta (batch 1): 0.229788
#> betas:    1.3151  -0.2309   0.3029   0.8384   0.5088   0.8119   0.3173  -0.1533   0.4592  -0.0091  -0.3327   0.1372
#> 
#> sigma: 0.147346
#> lambda (gene 1): 3.51505
#> epsilon (gene 1): -0.243696
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 655.128
#> delta (gene 1): 0.957524
#> phi (cell 1): 0.679813
#> s (cell 1): 0.403546
#> nu (cell 1): 0.409406
#> theta (batch 1): 0.268415
#> betas:    0.6891  -0.1027   1.1842   0.0766   0.8596   0.9329   0.4123   0.1493  -0.0202   0.2121  -0.3108  -0.3890
#> 
#> sigma: 0.19818
#> lambda (gene 1): 0.751749
#> epsilon (gene 1): 0.0552077
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 750.189
#> delta (gene 1): 0.465231
#> phi (cell 1): 0.691885
#> s (cell 1): 0.466868
#> nu (cell 1): 0.406835
#> theta (batch 1): 0.166709
#> betas:    0.6429  -0.0047   0.9277   0.6096   0.6314   0.8730   0.0406   0.0999  -0.2052  -0.2038  -0.4613  -0.4554
#> 
#> sigma: 0.125668
#> lambda (gene 1): 0.694529
#> epsilon (gene 1): -0.581333
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 933.067
#> delta (gene 1): 0.654091
#> phi (cell 1): 0.651739
#> s (cell 1): 0.522952
#> nu (cell 1): 0.405427
#> theta (batch 1): 0.142302
#> betas:    0.9699  -0.2103   1.3670   0.3964   0.7386   0.4740   0.5716  -0.0188   0.4461  -0.3930   0.4746  -0.2926
#> 
#> sigma: 0.159691
#> lambda (gene 1): 0.38688
#> epsilon (gene 1): -0.186834
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 623.928
#> delta (gene 1): 0.759287
#> phi (cell 1): 0.621703
#> s (cell 1): 0.325195
#> nu (cell 1): 0.413139
#> theta (batch 1): 0.073554
#> betas:    0.8100   0.0152   0.7915   1.0229   0.3474   0.5197   0.3150  -0.1048  -0.1225  -0.4153  -0.8627  -0.3704
#> 
#> sigma: 0.136479
#> lambda (gene 1): 0.627302
#> epsilon (gene 1): 0.0818124
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 704.015
#> delta (gene 1): 0.612437
#> phi (cell 1): 0.663294
#> s (cell 1): 0.48163
#> nu (cell 1): 0.414047
#> theta (batch 1): 0.024214
#> betas:    1.0593  -0.0841   0.7137   0.6350   0.8576   0.3290   0.2913  -0.0255  -0.0343  -0.2987  -0.2651  -0.7640
#> 
#> sigma: 0.171443
#> lambda (gene 1): 1.00944
#> epsilon (gene 1): -0.251119
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 540.3
#> delta (gene 1): 0.572138
#> phi (cell 1): 0.72685
#> s (cell 1): 0.330647
#> nu (cell 1): 0.438244
#> theta (batch 1): 0.0177877
#> betas:    1.3264  -0.0884   0.6879   0.5139   0.0775   0.3515   0.5247  -0.4876  -0.0075  -0.1004  -0.7895  -0.4636
#> 
#> sigma: 0.179416
#> lambda (gene 1): 0.871929
#> epsilon (gene 1): -0.370499
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 843.044
#> delta (gene 1): 0.737537
#> phi (cell 1): 0.721312
#> s (cell 1): 0.417875
#> nu (cell 1): 0.425244
#> theta (batch 1): 0.0315364
#> betas:    1.0248  -0.0705   0.5951   1.2219   0.1830   0.5175   0.2648   0.2162  -0.2675  -0.4884  -0.3730  -0.4075
#> 
#> sigma: 0.166869
#> lambda (gene 1): 1.3274
#> epsilon (gene 1): 0.146884
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 675.739
#> delta (gene 1): 0.87517
#> phi (cell 1): 0.649152
#> s (cell 1): 0.578951
#> nu (cell 1): 0.427094
#> theta (batch 1): 0.162876
#> betas:    0.9758  -0.1420   0.6773   1.0811   0.5814   0.2477   0.6915  -0.3476   0.6858  -0.1194  -0.5802   0.2643
#> 
#> sigma: 0.196483
#> lambda (gene 1): 1.36652
#> epsilon (gene 1): -0.0757794
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 735.502
#> delta (gene 1): 0.756056
#> phi (cell 1): 0.708538
#> s (cell 1): 0.419288
#> nu (cell 1): 0.420385
#> theta (batch 1): 0.105953
#> betas:    0.8571  -0.1257   0.4512   0.6793   0.6980   0.4664   0.4429   0.4122  -0.3655   0.1194   0.1064  -0.4714
#> 
#> sigma: 0.156516
#> lambda (gene 1): 1.31521
#> epsilon (gene 1): -0.22919
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 519.132
#> delta (gene 1): 1.12226
#> phi (cell 1): 0.651102
#> s (cell 1): 0.544342
#> nu (cell 1): 0.424576
#> theta (batch 1): 0.0448283
#> betas:    0.7382   0.0021   1.2714   0.3009   0.6658   0.3027   0.6195  -0.3411   0.3218  -0.3438  -0.8011  -0.2293
#> 
#> sigma: 0.156176
#> lambda (gene 1): 1.29491
#> epsilon (gene 1): 0.224052
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 542.38
#> delta (gene 1): 0.615602
#> phi (cell 1): 0.552977
#> s (cell 1): 0.173797
#> nu (cell 1): 0.427228
#> theta (batch 1): 0.543285
#> betas:    0.9882  -0.0944   0.4710   0.5824   0.5557   0.7581   0.3116   0.1806  -0.1763  -0.0792  -0.2706  -0.5107
#> 
#> sigma: 0.129769
#> lambda (gene 1): 0.820365
#> epsilon (gene 1): -0.347522
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 654.191
#> delta (gene 1): 0.713211
#> phi (cell 1): 0.6494
#> s (cell 1): 1.2313
#> nu (cell 1): 0.431185
#> theta (batch 1): 0.247327
#> betas:    0.7727  -0.0816   0.3238   0.8358   0.8149   0.9958  -0.1139   0.2854   0.3830  -0.5167  -0.1868  -0.3673
#> 
#> sigma: 0.154399
#> lambda (gene 1): 1.1219
#> epsilon (gene 1): -0.111324
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 622.489
#> delta (gene 1): 0.624246
#> phi (cell 1): 0.707974
#> s (cell 1): 0.582799
#> nu (cell 1): 0.413236
#> theta (batch 1): 0.460936
#> betas:    0.2680  -0.0574   0.7174   0.7578   1.0779   0.5329   0.4995   0.1703   0.4570  -0.2158  -0.3824   0.0997
#> 
#> sigma: 0.138703
#> lambda (gene 1): 2.0646
#> epsilon (gene 1): -0.230863
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 667.133
#> delta (gene 1): 0.575407
#> phi (cell 1): 0.734488
#> s (cell 1): 0.422351
#> nu (cell 1): 0.430807
#> theta (batch 1): 0.316906
#> betas:    0.7864   0.0574   0.8410   0.7200   0.8715  -0.0558   0.3692  -0.3116  -0.2280  -0.3815  -0.6021  -1.0359
#> 
#> sigma: 0.181341
#> lambda (gene 1): 2.53097
#> epsilon (gene 1): -0.378993
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.246
#> delta (gene 1): 0.712637
#> phi (cell 1): 0.78567
#> s (cell 1): 0.353413
#> nu (cell 1): 0.435543
#> theta (batch 1): 0.228821
#> betas:    0.6778  -0.0465   0.3317   0.9880   0.6589   0.5456   0.5389   0.0866  -0.0147  -0.3067  -0.1460  -0.5753
#> 
#> sigma: 0.121905
#> lambda (gene 1): 1.38661
#> epsilon (gene 1): -0.151717
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.465
#> delta (gene 1): 0.569111
#> phi (cell 1): 0.663532
#> s (cell 1): 0.546418
#> nu (cell 1): 0.415108
#> theta (batch 1): 0.112224
#> betas:    0.9100  -0.0193   1.1331  -0.2316   1.0656   0.5645  -0.2572   0.4340  -0.2027  -1.0127   0.0330  -0.6815
#> 
#> sigma: 0.187146
#> lambda (gene 1): 0.478577
#> epsilon (gene 1): -0.116352
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 605.049
#> delta (gene 1): 0.980757
#> phi (cell 1): 0.680836
#> s (cell 1): 0.351992
#> nu (cell 1): 0.419346
#> theta (batch 1): 0.0491752
#> betas:    0.8095  -0.1193   0.4466   0.8176   0.6062   0.7097   0.4419   0.2491  -0.0039  -0.3624  -0.0553  -0.2478
#> 
#> sigma: 0.140488
#> lambda (gene 1): 0.439329
#> epsilon (gene 1): 0.372915
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 573.939
#> delta (gene 1): 1.00913
#> phi (cell 1): 0.85944
#> s (cell 1): 0.466545
#> nu (cell 1): 0.413322
#> theta (batch 1): 0.0205696
#> betas:    0.6790  -0.1759   0.6931   1.0477   0.6860   0.8664   0.6177  -0.0984   0.3470   0.0251   0.2148  -0.0330
#> 
#> sigma: 0.132091
#> lambda (gene 1): 0.90818
#> epsilon (gene 1): 0.0933348
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.341
#> delta (gene 1): 0.646781
#> phi (cell 1): 0.804052
#> s (cell 1): 0.424334
#> nu (cell 1): 0.415304
#> theta (batch 1): 0.040783
#> betas:    0.1125  -0.0278   1.3598   1.0499   0.9249   0.3810   0.4563   0.2159   0.3360  -0.3233  -0.2189  -0.0284
#> 
#> sigma: 0.144701
#> lambda (gene 1): 0.46462
#> epsilon (gene 1): -0.146844
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 611.219
#> delta (gene 1): 0.701986
#> phi (cell 1): 0.726078
#> s (cell 1): 0.459638
#> nu (cell 1): 0.417217
#> theta (batch 1): 0.0134739
#> betas:    0.3677   0.0421   1.0254   0.8296   0.9327   0.4157   0.0856   0.3153  -0.1116  -0.5199  -0.1381  -0.9148
#> 
#> sigma: 0.172586
#> lambda (gene 1): 1.19661
#> epsilon (gene 1): -0.110201
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 773.554
#> delta (gene 1): 0.474654
#> phi (cell 1): 0.751152
#> s (cell 1): 0.435815
#> nu (cell 1): 0.427903
#> theta (batch 1): 0.00899034
#> betas:    1.4783  -0.3578   0.4726   0.5772   0.8893   0.3768   0.2525   0.7235   0.1368   0.4236  -0.2613   0.8419
#> 
#> sigma: 0.176208
#> lambda (gene 1): 0.25754
#> epsilon (gene 1): -0.561918
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 681.304
#> delta (gene 1): 0.794059
#> phi (cell 1): 0.692305
#> s (cell 1): 0.430264
#> nu (cell 1): 0.424637
#> theta (batch 1): 0.018722
#> betas:    0.5466   0.0980   0.5959  -0.1685   1.6167   0.0292   0.4094  -0.2973  -0.0353  -0.7348  -0.3804  -1.0105
#> 
#> sigma: 0.184296
#> lambda (gene 1): 3.09806
#> epsilon (gene 1): -0.041602
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 707.02
#> delta (gene 1): 0.946245
#> phi (cell 1): 0.69751
#> s (cell 1): 0.585266
#> nu (cell 1): 0.414185
#> theta (batch 1): 0.0301935
#> betas:    0.3147  -0.0823   1.1800   0.6025   0.5698   1.1244   0.3541   0.2857   0.4854  -0.2209  -0.3562   0.1634
#> 
#> sigma: 0.159137
#> lambda (gene 1): 0.927656
#> epsilon (gene 1): 0.231016
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 767.086
#> delta (gene 1): 0.955166
#> phi (cell 1): 0.73068
#> s (cell 1): 0.450364
#> nu (cell 1): 0.435965
#> theta (batch 1): 0.0797948
#> betas:    1.2743  -0.1568   0.7369   0.5884   0.3968   0.8630  -0.0236   0.1722   0.2905  -0.7557   0.1865  -0.3873
#> 
#> sigma: 0.164926
#> lambda (gene 1): 1.19902
#> epsilon (gene 1): 0.219922
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 770.982
#> delta (gene 1): 0.566931
#> phi (cell 1): 0.779359
#> s (cell 1): 0.3297
#> nu (cell 1): 0.416212
#> theta (batch 1): 0.0517453
#> betas:    0.4922   0.1452   1.1731   0.6493  -0.0222   0.6292   0.4762  -0.0756  -0.1782  -1.0040  -0.3391  -1.1414
#> 
#> sigma: 0.190226
#> lambda (gene 1): 1.01452
#> epsilon (gene 1): -0.321164
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 612.989
#> delta (gene 1): 0.642681
#> phi (cell 1): 0.890432
#> s (cell 1): 0.728485
#> nu (cell 1): 0.421076
#> theta (batch 1): 0.146937
#> betas:   -0.2076  -0.0618   0.5368   1.4467   0.9830   0.6338   0.8149   0.3357   0.5993  -0.1014   0.3651  -0.6334
#> 
#> sigma: 0.184031
#> lambda (gene 1): 1.52228
#> epsilon (gene 1): -0.305393
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 742.416
#> delta (gene 1): 0.544518
#> phi (cell 1): 0.743142
#> s (cell 1): 0.503569
#> nu (cell 1): 0.423008
#> theta (batch 1): 0.0299467
#> betas:    0.4211   0.0481   1.0287   1.0985   0.3783   0.4157   0.2646   0.7215  -0.7245  -0.3196  -0.3282  -0.2862
#> 
#> sigma: 0.188299
#> lambda (gene 1): 0.682462
#> epsilon (gene 1): -0.389599
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 642.052
#> delta (gene 1): 0.802123
#> phi (cell 1): 0.794029
#> s (cell 1): 0.478099
#> nu (cell 1): 0.435557
#> theta (batch 1): 0.11731
#> betas:    0.9365  -0.1704   0.0856   1.0625   0.6269   0.3000   1.0439  -0.1644   0.0244   0.0138   0.2020  -0.3445
#> 
#> sigma: 0.115568
#> lambda (gene 1): 0.419096
#> epsilon (gene 1): -0.108032
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 706.008
#> delta (gene 1): 0.708699
#> phi (cell 1): 0.885749
#> s (cell 1): 0.435446
#> nu (cell 1): 0.40605
#> theta (batch 1): 0.0175975
#> betas:    0.9108  -0.0572   0.9696   0.4181   0.5407   0.8039  -0.0255   0.1814  -0.1214  -0.2490  -0.3142  -0.5157
#> 
#> sigma: 0.129938
#> lambda (gene 1): 0.915277
#> epsilon (gene 1): -0.197912
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 700.245
#> delta (gene 1): 0.895379
#> phi (cell 1): 0.849639
#> s (cell 1): 0.352141
#> nu (cell 1): 0.420775
#> theta (batch 1): 0.02274
#> betas:    1.4235  -0.2190  -0.1469   0.5273   1.2219  -0.0823   0.6612   0.0331  -0.3729   0.1835  -0.2622  -0.1111
#> 
#> sigma: 0.161441
#> lambda (gene 1): 0.808317
#> epsilon (gene 1): 0.185566
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 11980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 750.705
#> delta (gene 1): 0.697322
#> phi (cell 1): 0.818165
#> s (cell 1): 0.477976
#> nu (cell 1): 0.425486
#> theta (batch 1): 0.0107537
#> betas:    0.8785  -0.1182   0.9229   0.3107   1.1086   0.4135   0.2460   0.2956   0.2790  -0.2510  -0.2506  -0.1226
#> 
#> sigma: 0.12295
#> lambda (gene 1): 0.799655
#> epsilon (gene 1): -0.228628
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 583.919
#> delta (gene 1): 0.707409
#> phi (cell 1): 0.77054
#> s (cell 1): 0.37752
#> nu (cell 1): 0.431503
#> theta (batch 1): 0.0818308
#> betas:    1.2102  -0.0996   1.1809   0.2526   0.2711   1.1408  -0.8327   0.6171  -0.2498  -0.0656  -0.6457   0.0394
#> 
#> sigma: 0.228903
#> lambda (gene 1): 2.15271
#> epsilon (gene 1): -0.290458
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 522.656
#> delta (gene 1): 1.07199
#> phi (cell 1): 0.685846
#> s (cell 1): 0.377009
#> nu (cell 1): 0.411799
#> theta (batch 1): 0.0760774
#> betas:    0.6118  -0.0760   0.7279   1.2355   0.3935   0.3915   0.5599   0.6617  -0.6197   0.3203  -0.2777  -0.2148
#> 
#> sigma: 0.162838
#> lambda (gene 1): 1.20826
#> epsilon (gene 1): 0.161336
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 519.15
#> delta (gene 1): 0.997136
#> phi (cell 1): 0.844322
#> s (cell 1): 1.81607
#> nu (cell 1): 0.420055
#> theta (batch 1): 0.25793
#> betas:    0.7933  -0.0170   0.2893   0.9340   1.0651   0.2405  -0.1750   0.5934  -0.0113  -0.6015  -0.1556  -0.7757
#> 
#> sigma: 0.115508
#> lambda (gene 1): 1.89999
#> epsilon (gene 1): 0.147966
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 697.949
#> delta (gene 1): 0.552708
#> phi (cell 1): 0.917822
#> s (cell 1): 1.22981
#> nu (cell 1): 0.428666
#> theta (batch 1): 0.219535
#> betas:    1.0198  -0.0998   0.1687   0.9752   0.7674   0.3258   0.3123   0.1677   0.1740  -0.0935  -0.8724   0.1331
#> 
#> sigma: 0.171615
#> lambda (gene 1): 1.37952
#> epsilon (gene 1): -0.38946
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 632.597
#> delta (gene 1): 0.453286
#> phi (cell 1): 0.981992
#> s (cell 1): 0.548924
#> nu (cell 1): 0.405606
#> theta (batch 1): 0.376818
#> betas:    0.9608  -0.1921   0.9175   0.5586   1.0599   0.4791   0.2757   0.4280  -0.0231  -0.0516   0.1018  -0.0408
#> 
#> sigma: 0.206612
#> lambda (gene 1): 1.38663
#> epsilon (gene 1): -0.574535
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 654.605
#> delta (gene 1): 0.528841
#> phi (cell 1): 0.849301
#> s (cell 1): 0.647498
#> nu (cell 1): 0.432904
#> theta (batch 1): 0.179384
#> betas:    0.9211  -0.0178   1.0455   0.5993   0.5192   0.6569  -0.2772   0.6224  -0.7167  -0.1378  -0.6915  -0.3714
#> 
#> sigma: 0.165469
#> lambda (gene 1): 0.743895
#> epsilon (gene 1): -0.325817
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 612.068
#> delta (gene 1): 0.507046
#> phi (cell 1): 0.642409
#> s (cell 1): 0.217719
#> nu (cell 1): 0.414754
#> theta (batch 1): 0.323812
#> betas:    0.3554   0.0676   0.6758   0.7101   0.6519   0.6907   0.3037  -0.3047   0.0900  -0.0529  -1.1192  -0.2902
#> 
#> sigma: 0.168518
#> lambda (gene 1): 1.75935
#> epsilon (gene 1): -0.500118
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 535.323
#> delta (gene 1): 0.938219
#> phi (cell 1): 0.676142
#> s (cell 1): 1.8735
#> nu (cell 1): 0.43186
#> theta (batch 1): 0.557344
#> betas:    1.1284  -0.1124   0.8570   0.9122   0.4450   0.5858  -0.3187   0.8751  -0.5209  -0.0229  -0.3127  -0.3756
#> 
#> sigma: 0.227136
#> lambda (gene 1): 2.07239
#> epsilon (gene 1): 0.0891232
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 688.182
#> delta (gene 1): 0.669444
#> phi (cell 1): 0.641841
#> s (cell 1): 0.375623
#> nu (cell 1): 0.432429
#> theta (batch 1): 0.560638
#> betas:    0.9101  -0.0773   0.3131   0.4310   0.9025   0.4162  -0.0377   0.5765  -0.3383  -0.0864  -0.4421  -0.4066
#> 
#> sigma: 0.209661
#> lambda (gene 1): 1.22238
#> epsilon (gene 1): -0.15752
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 703.351
#> delta (gene 1): 0.778171
#> phi (cell 1): 0.653481
#> s (cell 1): 0.210134
#> nu (cell 1): 0.417671
#> theta (batch 1): 0.291285
#> betas:    0.3871   0.0804   0.3629   0.7478   1.1813   0.1444   0.6912  -0.1145  -0.6109   0.0661  -0.6452  -0.8747
#> 
#> sigma: 0.187741
#> lambda (gene 1): 0.475384
#> epsilon (gene 1): -0.0953284
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 780.096
#> delta (gene 1): 0.934188
#> phi (cell 1): 0.752985
#> s (cell 1): 0.746796
#> nu (cell 1): 0.43005
#> theta (batch 1): 0.173619
#> betas:    1.0827  -0.1330   1.0047   0.5277   0.4454   0.9918   0.1106  -0.1773   0.6518  -0.6927  -0.0609   0.0694
#> 
#> sigma: 0.15757
#> lambda (gene 1): 0.955837
#> epsilon (gene 1): 0.0576332
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 790.278
#> delta (gene 1): 0.637632
#> phi (cell 1): 0.738447
#> s (cell 1): 0.44814
#> nu (cell 1): 0.440401
#> theta (batch 1): 0.115574
#> betas:    0.6112  -0.0444   0.6611   0.8620   0.6526   0.4786   0.4601   0.2636  -0.2766  -0.1390  -0.2786  -0.5537
#> 
#> sigma: 0.108658
#> lambda (gene 1): 0.920805
#> epsilon (gene 1): -0.147094
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 801.532
#> delta (gene 1): 0.637825
#> phi (cell 1): 0.758755
#> s (cell 1): 0.467032
#> nu (cell 1): 0.414361
#> theta (batch 1): 0.201131
#> betas:    0.9813  -0.2372   0.6756   0.4870   0.8499   0.2160   0.7840   0.3957   0.1091   0.2709   0.1607  -0.0627
#> 
#> sigma: 0.156388
#> lambda (gene 1): 0.885088
#> epsilon (gene 1): -0.408136
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 924.77
#> delta (gene 1): 0.656932
#> phi (cell 1): 0.815179
#> s (cell 1): 0.580947
#> nu (cell 1): 0.441371
#> theta (batch 1): 0.0951833
#> betas:    0.8001  -0.0275   0.3799   0.7488   0.2947   0.5945   0.2197   0.0004   0.2664  -0.6730  -0.4018  -0.3340
#> 
#> sigma: 0.13664
#> lambda (gene 1): 0.551371
#> epsilon (gene 1): -0.14565
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 557.537
#> delta (gene 1): 1.00441
#> phi (cell 1): 0.758368
#> s (cell 1): 0.591222
#> nu (cell 1): 0.424524
#> theta (batch 1): 0.236262
#> betas:    1.1351  -0.0416   0.2142   0.8716   0.6248   0.0192   0.5364  -0.1422  -0.0124  -0.5158  -0.5295  -0.5622
#> 
#> sigma: 0.162498
#> lambda (gene 1): 0.849233
#> epsilon (gene 1): 0.247067
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 808.23
#> delta (gene 1): 0.609005
#> phi (cell 1): 0.703994
#> s (cell 1): 0.583983
#> nu (cell 1): 0.424073
#> theta (batch 1): 0.0260951
#> betas:    0.8535  -0.0398   0.7443   1.1514   0.3951   0.7027   0.3025  -0.0017   0.1790  -0.6026  -0.5106  -0.2680
#> 
#> sigma: 0.127939
#> lambda (gene 1): 1.34239
#> epsilon (gene 1): -0.145475
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 581.725
#> delta (gene 1): 0.745766
#> phi (cell 1): 0.689173
#> s (cell 1): 0.327837
#> nu (cell 1): 0.426613
#> theta (batch 1): 0.217164
#> betas:    0.8434  -0.0702   0.0074   0.9016   0.5941   0.4956   0.1040   0.2199  -0.1512  -0.0883  -0.6399   0.2461
#> 
#> sigma: 0.154753
#> lambda (gene 1): 1.78114
#> epsilon (gene 1): -0.119118
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 551.307
#> delta (gene 1): 0.642539
#> phi (cell 1): 0.573936
#> s (cell 1): 0.539151
#> nu (cell 1): 0.428737
#> theta (batch 1): 0.214082
#> betas:    1.0492  -0.0232   0.7205   0.8815  -0.0778   1.1546  -0.0381  -0.2665  -0.1210  -0.6042   0.0417  -1.0801
#> 
#> sigma: 0.224326
#> lambda (gene 1): 0.575327
#> epsilon (gene 1): -0.31984
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 966.903
#> delta (gene 1): 0.619286
#> phi (cell 1): 0.631649
#> s (cell 1): 0.378153
#> nu (cell 1): 0.408488
#> theta (batch 1): 0.331162
#> betas:    0.9148  -0.0531   0.0351   0.7952   0.4070   0.8226   0.0288  -0.0205   0.2332  -0.6409  -0.3686  -0.4860
#> 
#> sigma: 0.186734
#> lambda (gene 1): 2.17768
#> epsilon (gene 1): -0.130421
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 782.714
#> delta (gene 1): 0.607854
#> phi (cell 1): 0.675729
#> s (cell 1): 0.623455
#> nu (cell 1): 0.405714
#> theta (batch 1): 0.189108
#> betas:    1.1386  -0.0730   1.1786   0.6157   0.6239   0.0573   0.3896   0.0587   0.0198  -0.4123  -0.5655  -0.2973
#> 
#> sigma: 0.158231
#> lambda (gene 1): 0.519299
#> epsilon (gene 1): -0.261295
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 659.27
#> delta (gene 1): 0.669789
#> phi (cell 1): 0.722538
#> s (cell 1): 0.468577
#> nu (cell 1): 0.424574
#> theta (batch 1): 0.216178
#> betas:    1.0938  -0.0943   0.9186  -0.1502   1.0090   0.4884   0.0695  -0.0660   0.3082  -0.1779  -0.9178  -0.0008
#> 
#> sigma: 0.169576
#> lambda (gene 1): 1.11821
#> epsilon (gene 1): -0.179711
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 750.715
#> delta (gene 1): 0.784375
#> phi (cell 1): 0.692767
#> s (cell 1): 0.270706
#> nu (cell 1): 0.418603
#> theta (batch 1): 0.357471
#> betas:    0.3838   0.0617   0.6043   0.5073   1.3581   0.3085   0.0405   0.0707  -0.1874  -0.5807  -0.3448  -0.6038
#> 
#> sigma: 0.157228
#> lambda (gene 1): 0.196206
#> epsilon (gene 1): 0.084371
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 640.291
#> delta (gene 1): 0.808705
#> phi (cell 1): 0.753936
#> s (cell 1): 1.34173
#> nu (cell 1): 0.426191
#> theta (batch 1): 0.214337
#> betas:    0.8620   0.0131   0.4530   0.1215   0.6994   0.6024   0.2990  -0.1593  -0.4227   0.1534  -1.1591  -0.2196
#> 
#> sigma: 0.247648
#> lambda (gene 1): 1.30466
#> epsilon (gene 1): -0.0700375
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 740.284
#> delta (gene 1): 0.684684
#> phi (cell 1): 0.730614
#> s (cell 1): 0.387405
#> nu (cell 1): 0.429797
#> theta (batch 1): 0.539375
#> betas:    0.9595  -0.1915   1.0375   0.3501   0.4346   1.1783   0.1470   0.3571   0.1805   0.0289  -0.1829   0.0582
#> 
#> sigma: 0.16161
#> lambda (gene 1): 1.42472
#> epsilon (gene 1): -0.171327
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 669.41
#> delta (gene 1): 0.808158
#> phi (cell 1): 0.781468
#> s (cell 1): 0.638869
#> nu (cell 1): 0.433377
#> theta (batch 1): 0.298076
#> betas:    0.3089  -0.0322   1.0552   0.8148   0.9304   0.4722   0.5236  -0.3784   0.5023  -0.0122  -0.7456  -0.0584
#> 
#> sigma: 0.181996
#> lambda (gene 1): 0.301676
#> epsilon (gene 1): 0.0246416
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 717.459
#> delta (gene 1): 0.906661
#> phi (cell 1): 0.736004
#> s (cell 1): 0.413991
#> nu (cell 1): 0.422986
#> theta (batch 1): 0.177183
#> betas:    0.6449  -0.0079   0.6954   0.8120   0.6410   0.5890  -0.0869   0.3889   0.0022  -0.6837  -0.1193  -0.3115
#> 
#> sigma: 0.154922
#> lambda (gene 1): 0.644033
#> epsilon (gene 1): 0.0900901
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 655.784
#> delta (gene 1): 0.629456
#> phi (cell 1): 0.743357
#> s (cell 1): 0.332898
#> nu (cell 1): 0.41972
#> theta (batch 1): 0.11974
#> betas:    0.6168  -0.1070   0.1378   1.0399   0.9035   0.2153   0.9377   0.0295  -0.0246  -0.1925   0.1181  -0.2665
#> 
#> sigma: 0.162363
#> lambda (gene 1): 0.526962
#> epsilon (gene 1): -0.211656
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 857.967
#> delta (gene 1): 0.850978
#> phi (cell 1): 0.719607
#> s (cell 1): 0.423123
#> nu (cell 1): 0.427589
#> theta (batch 1): 0.0570766
#> betas:    0.6944  -0.0245   0.1767   1.0309   0.2260   0.7010   0.4544  -0.0022  -0.1846  -0.4644   0.1227  -0.6741
#> 
#> sigma: 0.164466
#> lambda (gene 1): 0.972503
#> epsilon (gene 1): -0.0250721
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 634.591
#> delta (gene 1): 0.875521
#> phi (cell 1): 0.786637
#> s (cell 1): 0.365697
#> nu (cell 1): 0.41603
#> theta (batch 1): 0.0616823
#> betas:    0.3352   0.0480   0.1054   1.2830   0.8661   0.2658   0.5797   0.0697  -0.1226  -0.1260  -0.5870  -0.5872
#> 
#> sigma: 0.164205
#> lambda (gene 1): 0.952038
#> epsilon (gene 1): 0.000750823
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 737.131
#> delta (gene 1): 0.960631
#> phi (cell 1): 0.775561
#> s (cell 1): 0.488561
#> nu (cell 1): 0.417237
#> theta (batch 1): 0.0317632
#> betas:    0.7163  -0.0377   0.6735   0.8736   0.1964   0.6263   0.5482  -0.2916   0.3178  -0.5666   0.0850  -0.8474
#> 
#> sigma: 0.151041
#> lambda (gene 1): 1.0733
#> epsilon (gene 1): 0.075401
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 806.776
#> delta (gene 1): 0.735775
#> phi (cell 1): 0.768988
#> s (cell 1): 0.469477
#> nu (cell 1): 0.441392
#> theta (batch 1): 0.00905975
#> betas:    0.9094  -0.2332   0.4675   0.9502   1.0094   0.8644  -0.2190   0.7115   0.3259  -0.3232   0.5237   0.0976
#> 
#> sigma: 0.144954
#> lambda (gene 1): 1.35222
#> epsilon (gene 1): -0.127088
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 634.793
#> delta (gene 1): 0.841297
#> phi (cell 1): 0.781228
#> s (cell 1): 0.399034
#> nu (cell 1): 0.415672
#> theta (batch 1): 0.0519762
#> betas:    0.7390  -0.0360   0.4715   1.0665   0.2808   0.9643  -0.2151   0.6374  -0.3126  -0.3422  -0.4818  -0.1675
#> 
#> sigma: 0.183548
#> lambda (gene 1): 0.978633
#> epsilon (gene 1): 0.157253
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 802.654
#> delta (gene 1): 0.811771
#> phi (cell 1): 0.766043
#> s (cell 1): 0.579991
#> nu (cell 1): 0.391285
#> theta (batch 1): 0.070455
#> betas:    0.2042   0.0329   1.3330   0.4124   0.9514   0.4714   0.7326  -0.2226   0.0891   0.0907  -0.8095  -0.3846
#> 
#> sigma: 0.148829
#> lambda (gene 1): 0.727423
#> epsilon (gene 1): -0.0294873
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 692.89
#> delta (gene 1): 0.434667
#> phi (cell 1): 0.605155
#> s (cell 1): 0.430321
#> nu (cell 1): 0.431356
#> theta (batch 1): 0.106321
#> betas:    1.1797  -0.1507   0.9060   0.6901   0.1373   1.0030   0.1502   0.0357  -0.2472   0.0273   0.1672  -0.6613
#> 
#> sigma: 0.168128
#> lambda (gene 1): 0.795999
#> epsilon (gene 1): -0.826916
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 831.93
#> delta (gene 1): 0.727832
#> phi (cell 1): 0.668356
#> s (cell 1): 0.520115
#> nu (cell 1): 0.401226
#> theta (batch 1): 0.116167
#> betas:    0.8574  -0.1872   0.6889   0.7661   0.9859   0.6170   0.4712   0.1688   0.1768   0.2943  -0.4605   0.4060
#> 
#> sigma: 0.156172
#> lambda (gene 1): 0.828288
#> epsilon (gene 1): -0.15803
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 961.23
#> delta (gene 1): 0.758873
#> phi (cell 1): 0.706308
#> s (cell 1): 0.395991
#> nu (cell 1): 0.423788
#> theta (batch 1): 0.289553
#> betas:    0.8775  -0.0903   0.4484   0.9340   0.1458   0.4883   0.7105   0.0019  -0.0784  -0.1822   0.0026  -0.6068
#> 
#> sigma: 0.142592
#> lambda (gene 1): 1.24844
#> epsilon (gene 1): -0.171342
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 671.737
#> delta (gene 1): 0.73957
#> phi (cell 1): 0.64029
#> s (cell 1): 0.439455
#> nu (cell 1): 0.411452
#> theta (batch 1): 0.395843
#> betas:    0.6872  -0.0004   0.1994   1.1206   0.1597   0.8585   0.2254  -0.1963   0.0645  -0.3463  -0.3454  -0.4916
#> 
#> sigma: 0.125534
#> lambda (gene 1): 0.54705
#> epsilon (gene 1): -0.245346
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 711.389
#> delta (gene 1): 0.730725
#> phi (cell 1): 0.71157
#> s (cell 1): 0.398041
#> nu (cell 1): 0.41151
#> theta (batch 1): 0.161726
#> betas:    0.8550  -0.1415   0.8193   1.2098   0.2706   0.4844   0.4120   0.2692  -0.0185  -0.0329  -0.3690   0.1146
#> 
#> sigma: 0.163179
#> lambda (gene 1): 3.18493
#> epsilon (gene 1): -0.00659207
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 852.359
#> delta (gene 1): 0.56199
#> phi (cell 1): 0.796734
#> s (cell 1): 0.440488
#> nu (cell 1): 0.413274
#> theta (batch 1): 0.14188
#> betas:    0.6740  -0.0486   0.3973   0.6366   1.2107   0.3555   0.4367   0.1253   0.0475  -0.2726  -0.0796  -0.6595
#> 
#> sigma: 0.160523
#> lambda (gene 1): 0.772445
#> epsilon (gene 1): -0.47559
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 734.908
#> delta (gene 1): 0.628211
#> phi (cell 1): 0.731591
#> s (cell 1): 0.558919
#> nu (cell 1): 0.423589
#> theta (batch 1): 0.208686
#> betas:    0.9483  -0.0649   0.8670   0.4102   1.0258   0.0598   0.2456   0.2605   0.0678  -0.3893  -0.2983  -0.3454
#> 
#> sigma: 0.141257
#> lambda (gene 1): 0.365954
#> epsilon (gene 1): -0.373101
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 838.21
#> delta (gene 1): 0.831462
#> phi (cell 1): 0.728204
#> s (cell 1): 0.290347
#> nu (cell 1): 0.420387
#> theta (batch 1): 0.462334
#> betas:    0.8259  -0.1383   0.1045   0.9065   0.6706   0.6415   0.6623   0.0143   0.0354  -0.0054   0.0534  -0.2890
#> 
#> sigma: 0.151124
#> lambda (gene 1): 0.908767
#> epsilon (gene 1): -0.0917415
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 864.87
#> delta (gene 1): 0.757223
#> phi (cell 1): 0.675915
#> s (cell 1): 0.412572
#> nu (cell 1): 0.43852
#> theta (batch 1): 0.234827
#> betas:    0.7379  -0.1103   0.5562   0.8036   1.0547   0.4063   0.7157   0.0707   0.0930  -0.2870  -0.0302  -0.1812
#> 
#> sigma: 0.162709
#> lambda (gene 1): 1.28633
#> epsilon (gene 1): -0.0213806
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 928.327
#> delta (gene 1): 0.754174
#> phi (cell 1): 0.807656
#> s (cell 1): 0.566225
#> nu (cell 1): 0.414845
#> theta (batch 1): 0.275965
#> betas:    0.7023  -0.0009   0.7138   0.5881   0.6821   0.3246   0.2710   0.0794  -0.1882  -0.3616  -0.2573  -0.6862
#> 
#> sigma: 0.164169
#> lambda (gene 1): 1.80856
#> epsilon (gene 1): -0.127246
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 531.063
#> delta (gene 1): 0.607496
#> phi (cell 1): 0.880657
#> s (cell 1): 0.892442
#> nu (cell 1): 0.416173
#> theta (batch 1): 0.22084
#> betas:    1.1972  -0.2197   0.7515   0.5550   0.7933   0.6127  -0.0001   0.4788  -0.0196   0.0808  -0.1064   0.1258
#> 
#> sigma: 0.145642
#> lambda (gene 1): 0.731592
#> epsilon (gene 1): -0.359785
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 924.586
#> delta (gene 1): 1.01519
#> phi (cell 1): 0.996701
#> s (cell 1): 2.71427
#> nu (cell 1): 0.42181
#> theta (batch 1): 0.35207
#> betas:    0.9536  -0.0827   0.3889   0.8156   1.0034  -0.1510   0.6022   0.3733  -0.6835   0.0512  -0.2190  -0.5239
#> 
#> sigma: 0.217541
#> lambda (gene 1): 0.528866
#> epsilon (gene 1): 0.224974
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 741.652
#> delta (gene 1): 1.03142
#> phi (cell 1): 0.964976
#> s (cell 1): 1.06894
#> nu (cell 1): 0.424797
#> theta (batch 1): 0.178378
#> betas:    4.7366e-01  -8.5544e-02   8.2991e-01   7.2151e-01   5.3524e-01   6.1118e-01   7.3761e-01  -3.6927e-02   9.8104e-02   2.1675e-01  -4.8821e-01  -5.4725e-05
#> 
#> sigma: 0.128444
#> lambda (gene 1): 0.963726
#> epsilon (gene 1): 0.20018
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 817.705
#> delta (gene 1): 0.651027
#> phi (cell 1): 0.903266
#> s (cell 1): 0.463902
#> nu (cell 1): 0.410925
#> theta (batch 1): 0.161819
#> betas:    0.6811  -0.0107   0.5125   0.9470   0.3749   0.5406   0.2967   0.2156  -0.2816  -0.5419  -0.0652  -0.4701
#> 
#> sigma: 0.143106
#> lambda (gene 1): 2.71241
#> epsilon (gene 1): -0.187284
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 616.73
#> delta (gene 1): 0.573468
#> phi (cell 1): 0.868529
#> s (cell 1): 0.644647
#> nu (cell 1): 0.414294
#> theta (batch 1): 0.407141
#> betas:    0.4280   0.0283   1.0808   1.0207   0.2280   0.5802   0.9083  -0.4664  -0.0891  -0.0889  -0.4370  -0.6040
#> 
#> sigma: 0.154967
#> lambda (gene 1): 1.11995
#> epsilon (gene 1): -0.451109
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 564.865
#> delta (gene 1): 0.729595
#> phi (cell 1): 0.88308
#> s (cell 1): 1.18692
#> nu (cell 1): 0.403838
#> theta (batch 1): 0.345269
#> betas:    0.9303  -0.0307   0.9051   0.3600   0.8640   0.3640   0.3100  -0.0720  -0.0346  -0.5924   0.0520  -0.8848
#> 
#> sigma: 0.173321
#> lambda (gene 1): 0.839591
#> epsilon (gene 1): -0.207133
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 12980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 740.11
#> delta (gene 1): 0.515074
#> phi (cell 1): 0.893912
#> s (cell 1): 0.459771
#> nu (cell 1): 0.435196
#> theta (batch 1): 0.0820746
#> betas:    0.6222  -0.0738   0.4604   0.8342   1.0396   0.2947   0.0316   0.6529   0.0529  -0.1531  -0.5920  -0.0123
#> 
#> sigma: 0.150525
#> lambda (gene 1): 0.728695
#> epsilon (gene 1): -0.361565
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 750.366
#> delta (gene 1): 0.562971
#> phi (cell 1): 1.00779
#> s (cell 1): 0.226523
#> nu (cell 1): 0.420141
#> theta (batch 1): 0.458687
#> betas:    0.9887  -0.1384   0.1131   1.0059   0.6647   0.2637   0.5590  -0.0362  -0.1317   0.2582  -0.5204  -0.0503
#> 
#> sigma: 0.157401
#> lambda (gene 1): 0.520162
#> epsilon (gene 1): -0.414481
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 708.844
#> delta (gene 1): 0.47782
#> phi (cell 1): 1.09081
#> s (cell 1): 1.12869
#> nu (cell 1): 0.432234
#> theta (batch 1): 0.359526
#> betas:    1.3323  -0.1865  -0.0144   1.1841   0.3632   0.6461  -0.0100   0.4543  -0.2977  -0.2259   0.1187  -0.1131
#> 
#> sigma: 0.196203
#> lambda (gene 1): 0.221365
#> epsilon (gene 1): -0.567021
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1102.56
#> delta (gene 1): 0.850583
#> phi (cell 1): 1.14639
#> s (cell 1): 0.6064
#> nu (cell 1): 0.4151
#> theta (batch 1): 0.246879
#> betas:    1.2034  -0.1466   0.8249   0.2730   0.5764   0.8315  -0.0042  -0.1092   0.5648  -0.4393  -0.3690  -0.0619
#> 
#> sigma: 0.16738
#> lambda (gene 1): 1.06499
#> epsilon (gene 1): 0.00136096
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 861.352
#> delta (gene 1): 0.427553
#> phi (cell 1): 0.953661
#> s (cell 1): 0.376107
#> nu (cell 1): 0.412042
#> theta (batch 1): 0.216497
#> betas:   -0.2986   0.1399   0.6506   0.9449   0.9813   0.4483   0.6167   0.3590  -0.1451  -0.3211  -0.2793  -0.9429
#> 
#> sigma: 0.170918
#> lambda (gene 1): 0.372514
#> epsilon (gene 1): -0.678317
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 746.788
#> delta (gene 1): 0.615276
#> phi (cell 1): 0.969648
#> s (cell 1): 0.449216
#> nu (cell 1): 0.425938
#> theta (batch 1): 0.107818
#> betas:    0.4883   0.0295   0.9603   0.1918   0.8221   0.6786   0.3535  -0.3321   0.3220  -0.6332  -0.4112  -0.4581
#> 
#> sigma: 0.137064
#> lambda (gene 1): 0.272485
#> epsilon (gene 1): -0.245996
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 543.988
#> delta (gene 1): 0.895292
#> phi (cell 1): 0.834018
#> s (cell 1): 0.371032
#> nu (cell 1): 0.430072
#> theta (batch 1): 0.0513434
#> betas:    0.6678  -0.1203   0.6077   1.1621   0.7695   0.5847   0.4481   0.4275   0.0342   0.0578  -0.6983   0.7132
#> 
#> sigma: 0.204444
#> lambda (gene 1): 1.19146
#> epsilon (gene 1): 0.159673
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 670.218
#> delta (gene 1): 0.766833
#> phi (cell 1): 0.949897
#> s (cell 1): 0.265049
#> nu (cell 1): 0.413371
#> theta (batch 1): 0.0873005
#> betas:    1.4671  -0.0168   0.2174   0.7421   0.5463   0.0539  -0.0024  -0.1286   0.0614  -0.9109  -0.3847  -0.9922
#> 
#> sigma: 0.169911
#> lambda (gene 1): 2.63638
#> epsilon (gene 1): -0.143998
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 695.808
#> delta (gene 1): 0.653761
#> phi (cell 1): 0.753127
#> s (cell 1): 0.518574
#> nu (cell 1): 0.42485
#> theta (batch 1): 0.30606
#> betas:    0.8900  -0.1116   0.9211   1.0054   0.2074   0.7113   0.7882  -0.2300   0.2145  -0.1232  -0.5448  -0.0432
#> 
#> sigma: 0.191813
#> lambda (gene 1): 0.776343
#> epsilon (gene 1): -0.141611
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 609.427
#> delta (gene 1): 0.63796
#> phi (cell 1): 0.814033
#> s (cell 1): 1.09272
#> nu (cell 1): 0.422047
#> theta (batch 1): 0.506605
#> betas:    0.6620  -0.1269   0.3179   0.8141   0.7278   0.6849   0.3636   0.3677   0.0174  -0.0558  -0.0277  -0.0277
#> 
#> sigma: 0.198178
#> lambda (gene 1): 1.62474
#> epsilon (gene 1): -0.28571
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 827.948
#> delta (gene 1): 0.636007
#> phi (cell 1): 0.775895
#> s (cell 1): 0.483536
#> nu (cell 1): 0.428354
#> theta (batch 1): 0.169523
#> betas:    0.7389  -0.1752   1.2887   0.5126   0.6449   0.8919   0.2351   0.3655   0.4716  -0.2093  -0.1820   0.1354
#> 
#> sigma: 0.165959
#> lambda (gene 1): 0.513509
#> epsilon (gene 1): -0.111895
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 843.206
#> delta (gene 1): 0.716143
#> phi (cell 1): 0.720897
#> s (cell 1): 0.952462
#> nu (cell 1): 0.417742
#> theta (batch 1): 0.572461
#> betas:    0.4461  -0.0665   0.9108   0.3904   0.9932   0.5717   0.5092   0.1179   0.0797  -0.2893   0.0730  -0.3425
#> 
#> sigma: 0.128648
#> lambda (gene 1): 0.652042
#> epsilon (gene 1): -0.101209
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 635.544
#> delta (gene 1): 0.761777
#> phi (cell 1): 0.680185
#> s (cell 1): 0.585161
#> nu (cell 1): 0.4248
#> theta (batch 1): 0.168366
#> betas:    1.2985  -0.1449   0.9485   0.3261   0.3568   0.5964   0.2453  -0.0972  -0.0121  -0.1988  -0.5127   0.0744
#> 
#> sigma: 0.16469
#> lambda (gene 1): 1.03174
#> epsilon (gene 1): -0.0211056
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 757.134
#> delta (gene 1): 0.56956
#> phi (cell 1): 0.777784
#> s (cell 1): 0.469318
#> nu (cell 1): 0.440976
#> theta (batch 1): 0.0640216
#> betas:    0.8784  -0.1451   0.9631   0.7542   0.8082   0.1734  -0.0199   0.8859   0.1053  -0.5648  -0.3002   0.4740
#> 
#> sigma: 0.166203
#> lambda (gene 1): 0.705709
#> epsilon (gene 1): -0.078605
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 625.808
#> delta (gene 1): 0.460643
#> phi (cell 1): 0.758345
#> s (cell 1): 0.493747
#> nu (cell 1): 0.427565
#> theta (batch 1): 0.128194
#> betas:    0.8538  -0.1175   0.2772   0.7472   0.6416   0.6246   0.5234   0.1386  -0.4326   0.4053  -0.4223  -0.1622
#> 
#> sigma: 0.14439
#> lambda (gene 1): 0.76187
#> epsilon (gene 1): -0.649888
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 737.11
#> delta (gene 1): 0.64543
#> phi (cell 1): 0.73114
#> s (cell 1): 0.545247
#> nu (cell 1): 0.423085
#> theta (batch 1): 0.181794
#> betas:    1.1079  -0.0996   0.3449   0.9694   0.5227   0.3904   0.1662   0.3395  -0.0061  -0.4063  -0.4101  -0.2295
#> 
#> sigma: 0.152384
#> lambda (gene 1): 1.14914
#> epsilon (gene 1): -0.173392
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 881.792
#> delta (gene 1): 0.691331
#> phi (cell 1): 0.722616
#> s (cell 1): 0.490398
#> nu (cell 1): 0.420816
#> theta (batch 1): 0.202593
#> betas:    0.9850  -0.1312   0.2510   0.7596   0.6265   0.5587   0.4482  -0.0789   0.0070   0.2195  -0.6564   0.2319
#> 
#> sigma: 0.202378
#> lambda (gene 1): 0.268186
#> epsilon (gene 1): -0.271654
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 582.516
#> delta (gene 1): 0.73552
#> phi (cell 1): 0.801317
#> s (cell 1): 0.33614
#> nu (cell 1): 0.416438
#> theta (batch 1): 0.126553
#> betas:    0.9911  -0.1528   1.5825   0.3366   0.2266   0.8728   0.3669   0.1864   0.2131  -0.1665  -0.0650  -0.3328
#> 
#> sigma: 0.169495
#> lambda (gene 1): 1.70578
#> epsilon (gene 1): -0.16152
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 618.224
#> delta (gene 1): 0.554777
#> phi (cell 1): 0.808343
#> s (cell 1): 0.341604
#> nu (cell 1): 0.415673
#> theta (batch 1): 0.0650464
#> betas:    0.4829  -0.1630   0.7277   1.1618   0.3989   1.0153   0.6867   0.3548  -0.0846   0.2652   0.0796  -0.1260
#> 
#> sigma: 0.168223
#> lambda (gene 1): 0.606307
#> epsilon (gene 1): -0.321758
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 707.054
#> delta (gene 1): 0.694657
#> phi (cell 1): 0.813593
#> s (cell 1): 0.816896
#> nu (cell 1): 0.41245
#> theta (batch 1): 0.260357
#> betas:    0.7345   0.0325   0.6488   0.7486   0.5298   0.7305   0.0568   0.1928  -0.1060  -0.4924  -0.4771  -1.0277
#> 
#> sigma: 0.189775
#> lambda (gene 1): 0.903208
#> epsilon (gene 1): -0.137247
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 675.578
#> delta (gene 1): 0.450418
#> phi (cell 1): 0.801422
#> s (cell 1): 0.206363
#> nu (cell 1): 0.41992
#> theta (batch 1): 0.263848
#> betas:    0.7269  -0.0302   0.4043   0.7541   0.7591   0.4179   0.0870   0.1632   0.1027  -0.4952  -0.0812  -0.4995
#> 
#> sigma: 0.21029
#> lambda (gene 1): 1.48229
#> epsilon (gene 1): -0.728423
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 725.06
#> delta (gene 1): 0.75147
#> phi (cell 1): 0.905148
#> s (cell 1): 0.186814
#> nu (cell 1): 0.415123
#> theta (batch 1): 0.133923
#> betas:    0.7568  -0.1583   0.9299   0.9862   0.6581   0.8671   0.1585   0.5605   0.0550  -0.0104  -0.1431  -0.0753
#> 
#> sigma: 0.16283
#> lambda (gene 1): 0.212713
#> epsilon (gene 1): -0.00971956
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1037.37
#> delta (gene 1): 0.970207
#> phi (cell 1): 0.824557
#> s (cell 1): 0.406424
#> nu (cell 1): 0.426395
#> theta (batch 1): 0.286098
#> betas:    0.6291  -0.0819   0.5252   0.9516   0.4326   1.0497   0.3103   0.1005  -0.0199  -0.0577  -0.1857  -0.2957
#> 
#> sigma: 0.136168
#> lambda (gene 1): 0.481638
#> epsilon (gene 1): 0.148731
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 721.415
#> delta (gene 1): 0.536641
#> phi (cell 1): 0.82961
#> s (cell 1): 0.309106
#> nu (cell 1): 0.420611
#> theta (batch 1): 0.164619
#> betas:    1.1365  -0.1186  -0.0552   1.1526   0.5769   0.2852   0.3889   0.1960  -0.7278   0.2905  -0.4125  -0.3666
#> 
#> sigma: 0.197318
#> lambda (gene 1): 0.736023
#> epsilon (gene 1): -0.420762
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 922.151
#> delta (gene 1): 0.693577
#> phi (cell 1): 0.822666
#> s (cell 1): 0.57982
#> nu (cell 1): 0.432365
#> theta (batch 1): 0.575755
#> betas:    0.9833   0.0847   1.4240   0.0567   1.1109  -0.0298   0.0187  -0.0983  -0.5970  -0.3616  -0.8875  -1.1814
#> 
#> sigma: 0.210742
#> lambda (gene 1): 0.754902
#> epsilon (gene 1): -0.124342
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 595.974
#> delta (gene 1): 0.708689
#> phi (cell 1): 0.761623
#> s (cell 1): 0.871504
#> nu (cell 1): 0.422179
#> theta (batch 1): 0.181732
#> betas:    0.0565   0.1032   1.0513   0.7921   0.9261   0.5009   0.5862  -0.2443   0.0335  -0.0445  -0.8968  -0.5157
#> 
#> sigma: 0.178305
#> lambda (gene 1): 1.59691
#> epsilon (gene 1): -0.200953
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 613.421
#> delta (gene 1): 0.678768
#> phi (cell 1): 0.823906
#> s (cell 1): 0.254519
#> nu (cell 1): 0.418249
#> theta (batch 1): 0.112769
#> betas:   -0.0155   0.0519   0.7654   0.8789   1.1038   0.4260   0.2783   0.4817  -0.0649  -0.0639  -0.4507  -0.3263
#> 
#> sigma: 0.149236
#> lambda (gene 1): 0.620899
#> epsilon (gene 1): -0.285614
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 706.781
#> delta (gene 1): 0.750134
#> phi (cell 1): 0.833865
#> s (cell 1): 0.249365
#> nu (cell 1): 0.422329
#> theta (batch 1): 0.276879
#> betas:    1.2449  -0.2228   0.3559   0.6300   0.9432   0.1450   0.7859  -0.1897   0.5817  -0.1448  -0.2079   0.2231
#> 
#> sigma: 0.119736
#> lambda (gene 1): 1.03795
#> epsilon (gene 1): -0.177203
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 754.933
#> delta (gene 1): 0.672269
#> phi (cell 1): 0.763749
#> s (cell 1): 0.486922
#> nu (cell 1): 0.439075
#> theta (batch 1): 0.51073
#> betas:    0.0190  -0.0210   1.2163   1.0512   0.9428   0.8933  -0.0996   0.9530   0.3605  -0.6131   0.1295  -0.1598
#> 
#> sigma: 0.164557
#> lambda (gene 1): 0.835455
#> epsilon (gene 1): -0.152434
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 551.102
#> delta (gene 1): 0.874987
#> phi (cell 1): 0.849583
#> s (cell 1): 0.272383
#> nu (cell 1): 0.429069
#> theta (batch 1): 0.362735
#> betas:    0.9310  -0.1516   0.5519   0.5832   0.8288   0.8196   0.1244   0.3593  -0.0360   0.1402  -0.3575   0.0826
#> 
#> sigma: 0.132166
#> lambda (gene 1): 0.95178
#> epsilon (gene 1): -0.00583272
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 750.362
#> delta (gene 1): 0.634359
#> phi (cell 1): 0.931171
#> s (cell 1): 0.691336
#> nu (cell 1): 0.433615
#> theta (batch 1): 0.202638
#> betas:    0.9486  -0.0839   0.3571   1.3428   0.4476   0.5821   0.1042   0.2181  -0.1700  -0.0065  -0.6008  -0.1360
#> 
#> sigma: 0.146753
#> lambda (gene 1): 0.24871
#> epsilon (gene 1): -0.282084
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 831.438
#> delta (gene 1): 0.924638
#> phi (cell 1): 0.998671
#> s (cell 1): 0.412421
#> nu (cell 1): 0.42002
#> theta (batch 1): 0.144297
#> betas:    1.4898  -0.2278   0.2798   0.4074   0.4429   0.6776   0.3348  -0.2707   0.5413  -0.4936  -0.0507  -0.0821
#> 
#> sigma: 0.128779
#> lambda (gene 1): 0.891216
#> epsilon (gene 1): 0.1944
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 696.325
#> delta (gene 1): 0.542817
#> phi (cell 1): 0.928259
#> s (cell 1): 0.943577
#> nu (cell 1): 0.414799
#> theta (batch 1): 0.324164
#> betas:    0.7053   0.0199   0.6122   0.1533   1.0622   0.4435  -0.0542   0.1844  -0.1912  -0.3581  -0.9096  -0.2114
#> 
#> sigma: 0.132829
#> lambda (gene 1): 1.40898
#> epsilon (gene 1): -0.251374
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 773.067
#> delta (gene 1): 0.827202
#> phi (cell 1): 0.838602
#> s (cell 1): 0.838088
#> nu (cell 1): 0.417439
#> theta (batch 1): 0.578923
#> betas:    0.9560  -0.1119   0.7469   0.6100   0.6001   0.4722   0.6010  -0.0273   0.2022   0.0584  -0.5006  -0.3771
#> 
#> sigma: 0.143283
#> lambda (gene 1): 2.74002
#> epsilon (gene 1): -0.114606
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 734.522
#> delta (gene 1): 0.476128
#> phi (cell 1): 0.853764
#> s (cell 1): 0.663103
#> nu (cell 1): 0.417825
#> theta (batch 1): 0.314689
#> betas:    1.3999  -0.1241  -0.1696   1.2866   0.1595   0.6329  -0.0449   0.0391   0.1550  -0.5577  -0.3002  -0.3740
#> 
#> sigma: 0.15082
#> lambda (gene 1): 1.25438
#> epsilon (gene 1): -0.531964
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 660.444
#> delta (gene 1): 0.745612
#> phi (cell 1): 0.88462
#> s (cell 1): 1.45962
#> nu (cell 1): 0.433049
#> theta (batch 1): 0.532018
#> betas:    0.1747   0.0690   1.1783   0.8085   0.8940   0.1367   0.5561   0.1027  -0.0512  -0.1681  -1.0658  -0.0080
#> 
#> sigma: 0.156684
#> lambda (gene 1): 2.13623
#> epsilon (gene 1): 0.0472182
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 626.214
#> delta (gene 1): 0.760283
#> phi (cell 1): 0.958779
#> s (cell 1): 0.811107
#> nu (cell 1): 0.41364
#> theta (batch 1): 0.657342
#> betas:    0.8225  -0.1575   1.2342   0.3800   0.9666   0.7277  -0.0072   0.5051   0.0138   0.2700  -0.4230   0.0939
#> 
#> sigma: 0.183303
#> lambda (gene 1): 1.09727
#> epsilon (gene 1): -0.139755
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 627.916
#> delta (gene 1): 0.660087
#> phi (cell 1): 0.806513
#> s (cell 1): 0.344892
#> nu (cell 1): 0.422095
#> theta (batch 1): 0.128192
#> betas:    0.6013  -0.0616   1.0350   0.5990   0.9387   0.5411   0.4663   0.0888   0.1874  -0.0869  -0.4253  -0.3386
#> 
#> sigma: 0.134216
#> lambda (gene 1): 0.191164
#> epsilon (gene 1): -0.255174
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 841.772
#> delta (gene 1): 0.842876
#> phi (cell 1): 0.680824
#> s (cell 1): 0.348513
#> nu (cell 1): 0.431822
#> theta (batch 1): 0.17513
#> betas:    0.7735  -0.1190   0.8985   0.5748   1.0206   0.3587   0.2339   0.4550   0.0837  -0.3152  -0.0994   0.0989
#> 
#> sigma: 0.144604
#> lambda (gene 1): 0.412621
#> epsilon (gene 1): 0.0487561
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 674.466
#> delta (gene 1): 0.755506
#> phi (cell 1): 0.662951
#> s (cell 1): 0.618269
#> nu (cell 1): 0.412647
#> theta (batch 1): 0.240615
#> betas:    0.4426  -0.0877   0.7341   1.2586   0.6722   0.5782   0.7697   0.2281   0.1408  -0.1177  -0.3803   0.1185
#> 
#> sigma: 0.204401
#> lambda (gene 1): 0.392327
#> epsilon (gene 1): 0.065667
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 789.863
#> delta (gene 1): 0.644589
#> phi (cell 1): 0.788286
#> s (cell 1): 0.602247
#> nu (cell 1): 0.418979
#> theta (batch 1): 0.61861
#> betas:    0.4954  -0.0253   0.7501   0.8664   0.7950   0.3793   0.4583   0.2765  -0.3393   0.4327  -1.0519   0.0156
#> 
#> sigma: 0.241413
#> lambda (gene 1): 1.13574
#> epsilon (gene 1): -0.254283
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 740.185
#> delta (gene 1): 0.662725
#> phi (cell 1): 0.761825
#> s (cell 1): 0.599603
#> nu (cell 1): 0.413567
#> theta (batch 1): 0.0690542
#> betas:    1.3522  -0.1728   0.3474   0.9711   0.2401   0.7427   0.1684   0.0101   0.2683  -0.2582  -0.4218   0.0686
#> 
#> sigma: 0.149089
#> lambda (gene 1): 0.563142
#> epsilon (gene 1): -0.238897
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 736.198
#> delta (gene 1): 0.594348
#> phi (cell 1): 0.734872
#> s (cell 1): 1.13643
#> nu (cell 1): 0.436725
#> theta (batch 1): 0.154083
#> betas:    0.4430   0.0420   0.6839   0.8068   0.4404   0.7260   0.1191   0.4029  -0.4407   0.0792  -0.7428  -0.6894
#> 
#> sigma: 0.153821
#> lambda (gene 1): 1.14264
#> epsilon (gene 1): -0.359394
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 869.762
#> delta (gene 1): 0.603059
#> phi (cell 1): 0.59769
#> s (cell 1): 0.676146
#> nu (cell 1): 0.419613
#> theta (batch 1): 0.0869196
#> betas:    0.3295   0.0241   0.3625   0.7515   0.9980   0.3444   0.4792   0.2518   0.0922  -0.8995   0.1393  -0.5784
#> 
#> sigma: 0.177471
#> lambda (gene 1): 1.42622
#> epsilon (gene 1): -0.1681
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 673.816
#> delta (gene 1): 0.777564
#> phi (cell 1): 0.665623
#> s (cell 1): 0.569654
#> nu (cell 1): 0.411953
#> theta (batch 1): 0.203685
#> betas:    1.0520  -0.1806   1.2784   0.5538  -0.0508   1.0192  -0.0314   0.3016   0.2346  -0.0256   0.0279  -0.2437
#> 
#> sigma: 0.165335
#> lambda (gene 1): 0.11159
#> epsilon (gene 1): -0.244172
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 651.493
#> delta (gene 1): 0.867907
#> phi (cell 1): 0.680542
#> s (cell 1): 0.359929
#> nu (cell 1): 0.444417
#> theta (batch 1): 0.343497
#> betas:    0.9133  -0.0916   1.0744   0.5321   0.6885  -0.0511   0.9006  -0.1252  -0.0717  -0.2993  -0.2245  -0.1893
#> 
#> sigma: 0.102562
#> lambda (gene 1): 2.49127
#> epsilon (gene 1): 0.121164
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 570.658
#> delta (gene 1): 0.61978
#> phi (cell 1): 0.712961
#> s (cell 1): 0.341525
#> nu (cell 1): 0.422532
#> theta (batch 1): 0.219614
#> betas:    1.1098  -0.0611   0.7231   0.7246   0.4772   0.2109   0.8811  -0.7968   0.3067  -0.3096  -0.8498   0.0421
#> 
#> sigma: 0.17053
#> lambda (gene 1): 0.918239
#> epsilon (gene 1): -0.299149
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 630.267
#> delta (gene 1): 0.719641
#> phi (cell 1): 0.703197
#> s (cell 1): 1.92084
#> nu (cell 1): 0.415039
#> theta (batch 1): 0.353557
#> betas:    0.5358   0.1251   0.6318   0.5940   0.7438   0.4855  -0.2646   0.0565  -0.1801  -0.4632  -0.8794  -0.6553
#> 
#> sigma: 0.133041
#> lambda (gene 1): 2.70441
#> epsilon (gene 1): -0.262113
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 755.062
#> delta (gene 1): 0.770184
#> phi (cell 1): 0.665475
#> s (cell 1): 0.629117
#> nu (cell 1): 0.435377
#> theta (batch 1): 0.156726
#> betas:    1.0333   0.0206   0.6791   1.0618  -0.1409   0.5662   0.3363  -0.3280  -0.3830  -0.2713  -0.4333  -1.2719
#> 
#> sigma: 0.188946
#> lambda (gene 1): 1.1975
#> epsilon (gene 1): -0.155014
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 629.894
#> delta (gene 1): 0.697179
#> phi (cell 1): 0.721848
#> s (cell 1): 0.506489
#> nu (cell 1): 0.418058
#> theta (batch 1): 0.188174
#> betas:    0.3485   0.0013   0.7443   0.9365   0.5085   0.6812   0.3182   0.2045  -0.0419  -0.3202  -0.1920  -0.3891
#> 
#> sigma: 0.108199
#> lambda (gene 1): 0.123228
#> epsilon (gene 1): -0.168349
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 13980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 775.148
#> delta (gene 1): 0.828303
#> phi (cell 1): 0.718766
#> s (cell 1): 0.233329
#> nu (cell 1): 0.412625
#> theta (batch 1): 0.440179
#> betas:    1.1921  -0.1331   0.5822   0.4471   0.6139   0.7204   0.2383  -0.2624   0.5652  -0.7266  -0.1652  -0.2893
#> 
#> sigma: 0.154436
#> lambda (gene 1): 0.512752
#> epsilon (gene 1): 0.117317
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 600.916
#> delta (gene 1): 0.635006
#> phi (cell 1): 0.772514
#> s (cell 1): 0.68698
#> nu (cell 1): 0.422131
#> theta (batch 1): 0.325844
#> betas:    0.9483  -0.0531   0.0325   1.1867   0.3891   0.5917   0.3188  -0.2403   0.5291  -0.4106  -0.5635  -0.4381
#> 
#> sigma: 0.163358
#> lambda (gene 1): 0.529335
#> epsilon (gene 1): -0.387756
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 674.196
#> delta (gene 1): 0.916512
#> phi (cell 1): 0.762453
#> s (cell 1): 0.66692
#> nu (cell 1): 0.444772
#> theta (batch 1): 0.128585
#> betas:    0.6330  -0.1351   0.9501   1.0406   0.1379   1.4048   0.3763   0.1438   0.3576  -0.4021   0.2648  -0.2602
#> 
#> sigma: 0.150598
#> lambda (gene 1): 1.37787
#> epsilon (gene 1): 0.173955
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 773.154
#> delta (gene 1): 0.526902
#> phi (cell 1): 0.73575
#> s (cell 1): 0.924986
#> nu (cell 1): 0.428918
#> theta (batch 1): 0.17739
#> betas:    0.4599  -0.0841   0.5845   1.1211   0.7123   0.5637   0.5951   0.1707   0.1551  -0.2976  -0.0860   0.0996
#> 
#> sigma: 0.176706
#> lambda (gene 1): 0.927864
#> epsilon (gene 1): -0.366973
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 605.828
#> delta (gene 1): 0.549625
#> phi (cell 1): 0.695072
#> s (cell 1): 0.794598
#> nu (cell 1): 0.427312
#> theta (batch 1): 0.142406
#> betas:    0.6296  -0.0447   0.6042   1.0057   0.0533   0.7760   0.2473   0.3446  -0.1750  -0.2653  -0.0127  -0.4571
#> 
#> sigma: 0.164117
#> lambda (gene 1): 0.507542
#> epsilon (gene 1): -0.493148
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 845.163
#> delta (gene 1): 0.490407
#> phi (cell 1): 0.683384
#> s (cell 1): 0.643243
#> nu (cell 1): 0.41389
#> theta (batch 1): 0.0823541
#> betas:    1.0493  -0.1826   0.4756   0.6081   0.7120   0.6615   0.2690  -0.1583   0.7211  -0.2554  -0.1028  -0.0899
#> 
#> sigma: 0.158489
#> lambda (gene 1): 0.873307
#> epsilon (gene 1): -0.633291
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 666.532
#> delta (gene 1): 0.742923
#> phi (cell 1): 0.72425
#> s (cell 1): 0.352714
#> nu (cell 1): 0.422282
#> theta (batch 1): 0.118876
#> betas:    0.4723  -0.0142   0.8127   0.8398   0.7228   0.5850   0.2291   0.1939   0.0478  -0.1011  -0.8544   0.4262
#> 
#> sigma: 0.156056
#> lambda (gene 1): 1.14294
#> epsilon (gene 1): -0.128702
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 553.629
#> delta (gene 1): 0.631083
#> phi (cell 1): 0.711017
#> s (cell 1): 0.475392
#> nu (cell 1): 0.418163
#> theta (batch 1): 0.205163
#> betas:   -0.1486   0.0477   0.4787   1.8876   0.8851   0.7213   0.4514   0.2654  -0.0774   0.1320  -0.7469   0.2696
#> 
#> sigma: 0.187674
#> lambda (gene 1): 0.766697
#> epsilon (gene 1): -0.293078
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 620.361
#> delta (gene 1): 0.573396
#> phi (cell 1): 0.772159
#> s (cell 1): 0.832146
#> nu (cell 1): 0.42719
#> theta (batch 1): 0.483845
#> betas:    1.4881  -0.2215   0.2349   0.4115   0.8128   0.6530  -0.0051   0.2511  -0.2434   0.2453  -0.2251  -0.1780
#> 
#> sigma: 0.129555
#> lambda (gene 1): 0.601325
#> epsilon (gene 1): -0.487867
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 720.552
#> delta (gene 1): 0.616529
#> phi (cell 1): 0.731877
#> s (cell 1): 1.12857
#> nu (cell 1): 0.428612
#> theta (batch 1): 0.313827
#> betas:    0.7501  -0.0759   0.4408   1.1429   0.5022   0.5653   0.4143  -0.0437   0.1190  -0.2379  -0.5976   0.1853
#> 
#> sigma: 0.156244
#> lambda (gene 1): 2.20802
#> epsilon (gene 1): -0.174018
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 686.238
#> delta (gene 1): 0.912102
#> phi (cell 1): 0.723864
#> s (cell 1): 0.337069
#> nu (cell 1): 0.426653
#> theta (batch 1): 0.108043
#> betas:    0.8771  -0.1435   1.3031   0.7504   0.6425   0.8391   0.1894   0.1508   0.5207  -0.4800  -0.1082   0.1243
#> 
#> sigma: 0.185916
#> lambda (gene 1): 2.16058
#> epsilon (gene 1): 0.131092
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 662.618
#> delta (gene 1): 0.559394
#> phi (cell 1): 0.820455
#> s (cell 1): 0.276985
#> nu (cell 1): 0.423615
#> theta (batch 1): 0.459875
#> betas:    0.7507  -0.0900   0.4210   1.0669   1.0618   0.2285   0.2199   0.2960  -0.0020  -0.0460  -0.4309  -0.0478
#> 
#> sigma: 0.165855
#> lambda (gene 1): 0.519552
#> epsilon (gene 1): -0.414079
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 691.274
#> delta (gene 1): 0.722229
#> phi (cell 1): 0.816749
#> s (cell 1): 0.37242
#> nu (cell 1): 0.422708
#> theta (batch 1): 0.175962
#> betas:    1.1455  -0.1831   0.4360   0.9568   0.4385   0.6576  -0.1040   0.4037   0.1726  -0.1748  -0.2654   0.1663
#> 
#> sigma: 0.143999
#> lambda (gene 1): 2.74445
#> epsilon (gene 1): -0.120511
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 641.456
#> delta (gene 1): 0.703144
#> phi (cell 1): 0.821443
#> s (cell 1): 0.305702
#> nu (cell 1): 0.437129
#> theta (batch 1): 0.472722
#> betas:    1.1174  -0.1196   1.2995   0.4649   0.1887   0.8203   0.0598   0.1561  -0.2459   0.0874  -0.7532  -0.1933
#> 
#> sigma: 0.159971
#> lambda (gene 1): 1.61158
#> epsilon (gene 1): -0.00775681
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 861.64
#> delta (gene 1): 0.700498
#> phi (cell 1): 0.818114
#> s (cell 1): 0.865755
#> nu (cell 1): 0.419851
#> theta (batch 1): 0.264569
#> betas:    0.9344  -0.1625   0.6359   0.7293   0.9325   0.7205   0.2310   0.0799  -0.0655   0.0446  -0.2451   0.4214
#> 
#> sigma: 0.180575
#> lambda (gene 1): 1.40997
#> epsilon (gene 1): -0.17393
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 762.974
#> delta (gene 1): 0.622806
#> phi (cell 1): 0.705886
#> s (cell 1): 0.588537
#> nu (cell 1): 0.424448
#> theta (batch 1): 0.188714
#> betas:    1.3024  -0.0729   0.4621   0.6025   0.9234  -0.1873   0.5673  -0.0760  -0.3996  -0.1668  -0.7883  -0.3356
#> 
#> sigma: 0.171838
#> lambda (gene 1): 1.39667
#> epsilon (gene 1): -0.169748
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 794.244
#> delta (gene 1): 0.785576
#> phi (cell 1): 0.675088
#> s (cell 1): 0.548812
#> nu (cell 1): 0.414448
#> theta (batch 1): 0.206465
#> betas:    0.8713  -0.0569   0.5828   0.5744   0.8389   0.6942  -0.2573   0.1846   0.3277  -0.6639  -0.1767  -0.4853
#> 
#> sigma: 0.144684
#> lambda (gene 1): 1.20939
#> epsilon (gene 1): -0.0418961
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 865.94
#> delta (gene 1): 0.880374
#> phi (cell 1): 0.769147
#> s (cell 1): 0.284128
#> nu (cell 1): 0.410462
#> theta (batch 1): 0.265295
#> betas:    0.4143  -0.0138   1.0200   1.1556   0.3763   0.6607   0.4636  -0.1087   0.1057  -0.1619  -0.6393   0.0509
#> 
#> sigma: 0.155037
#> lambda (gene 1): 1.2036
#> epsilon (gene 1): 0.120291
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 908.284
#> delta (gene 1): 0.869174
#> phi (cell 1): 0.755082
#> s (cell 1): 0.503844
#> nu (cell 1): 0.442454
#> theta (batch 1): 0.261158
#> betas:    0.0576   0.0507   0.2409   0.9240   1.0801   0.5262   0.3504   0.0248   0.2760  -0.6037   0.1137  -0.8104
#> 
#> sigma: 0.145438
#> lambda (gene 1): 1.16718
#> epsilon (gene 1): 0.0174343
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1079.6
#> delta (gene 1): 0.965986
#> phi (cell 1): 0.794565
#> s (cell 1): 0.403291
#> nu (cell 1): 0.411051
#> theta (batch 1): 0.190865
#> betas:    0.9964  -0.0659   0.3568   0.6044   0.7524   0.3808   0.3572  -0.0410  -0.0462  -0.2908  -0.6529  -0.0578
#> 
#> sigma: 0.140259
#> lambda (gene 1): 1.04787
#> epsilon (gene 1): 0.237264
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 802.099
#> delta (gene 1): 1.07804
#> phi (cell 1): 0.758897
#> s (cell 1): 0.401201
#> nu (cell 1): 0.422335
#> theta (batch 1): 0.0273498
#> betas:    0.9650  -0.1668   0.6094   0.4220   1.1302   0.5087   0.2637   0.4977  -0.2594   0.3437  -0.5734   0.2480
#> 
#> sigma: 0.127066
#> lambda (gene 1): 0.52232
#> epsilon (gene 1): 0.295504
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 754.324
#> delta (gene 1): 0.789597
#> phi (cell 1): 0.660142
#> s (cell 1): 0.403356
#> nu (cell 1): 0.417901
#> theta (batch 1): 0.00971846
#> betas:    1.0407  -0.0514   0.0700   0.7303   0.5244   0.4037   0.0609   0.1403  -0.0103  -0.3987  -0.3087  -0.4808
#> 
#> sigma: 0.21611
#> lambda (gene 1): 1.60613
#> epsilon (gene 1): -0.188879
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 928.314
#> delta (gene 1): 0.818199
#> phi (cell 1): 0.597364
#> s (cell 1): 0.521872
#> nu (cell 1): 0.418171
#> theta (batch 1): 0.105751
#> betas:    0.8247  -0.1696   0.2312   1.0334   0.4754   0.9669   0.2362   0.3090   0.4165  -0.2060   0.1133  -0.2096
#> 
#> sigma: 0.151717
#> lambda (gene 1): 0.93925
#> epsilon (gene 1): -0.0634177
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 532.64
#> delta (gene 1): 0.647369
#> phi (cell 1): 0.625548
#> s (cell 1): 0.619726
#> nu (cell 1): 0.438978
#> theta (batch 1): 0.220684
#> betas:    1.2784  -0.1892   0.7175   0.5452   0.8145   0.3543   0.2046   0.1698   0.1713  -0.1759  -0.3976   0.2785
#> 
#> sigma: 0.118702
#> lambda (gene 1): 1.01063
#> epsilon (gene 1): -0.212529
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 707.643
#> delta (gene 1): 0.626368
#> phi (cell 1): 0.624103
#> s (cell 1): 0.604956
#> nu (cell 1): 0.428393
#> theta (batch 1): 0.480569
#> betas:    0.6276  -0.1081   0.7485   0.6160   0.7627   0.5260   0.5727   0.2539  -0.0090  -0.0707   0.1165  -0.2569
#> 
#> sigma: 0.131789
#> lambda (gene 1): 0.409877
#> epsilon (gene 1): -0.388445
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 906.898
#> delta (gene 1): 0.526764
#> phi (cell 1): 0.677256
#> s (cell 1): 2.05764
#> nu (cell 1): 0.409831
#> theta (batch 1): 0.48167
#> betas:    0.7009  -0.0369   0.6556   0.6313   0.7138   0.5200   0.3902  -0.0799   0.0428  -0.0491  -0.0875  -0.9843
#> 
#> sigma: 0.160608
#> lambda (gene 1): 0.362755
#> epsilon (gene 1): -0.71912
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 658.232
#> delta (gene 1): 0.529
#> phi (cell 1): 0.720633
#> s (cell 1): 0.575847
#> nu (cell 1): 0.424871
#> theta (batch 1): 0.415654
#> betas:    1.0381  -0.1828   0.4341   0.4326   1.0097   0.8036   0.0660   0.3946  -0.1625   0.2213  -0.6082   0.3152
#> 
#> sigma: 0.168999
#> lambda (gene 1): 0.772173
#> epsilon (gene 1): -0.289451
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 731.091
#> delta (gene 1): 0.487379
#> phi (cell 1): 0.766848
#> s (cell 1): 0.675488
#> nu (cell 1): 0.431388
#> theta (batch 1): 0.203537
#> betas:    1.4444  -0.2450   0.8309   0.5716   0.3909   0.8124   0.2542  -0.0223   0.5024  -0.6174   0.5167  -0.2730
#> 
#> sigma: 0.18603
#> lambda (gene 1): 0.0837789
#> epsilon (gene 1): -0.565122
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 686.747
#> delta (gene 1): 0.697023
#> phi (cell 1): 0.914559
#> s (cell 1): 0.329673
#> nu (cell 1): 0.429238
#> theta (batch 1): 0.440471
#> betas:    1.3590  -0.0588   0.0315   0.7373   0.5965   0.4086  -0.2708   0.2784  -0.2003  -0.6064  -0.4474  -0.4122
#> 
#> sigma: 0.17321
#> lambda (gene 1): 0.297673
#> epsilon (gene 1): -0.184453
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 647.142
#> delta (gene 1): 0.859274
#> phi (cell 1): 0.925754
#> s (cell 1): 0.395363
#> nu (cell 1): 0.413979
#> theta (batch 1): 0.276132
#> betas:    0.6645  -0.0920   0.1679   1.4120   0.3843   0.7515   0.5433  -0.4071   0.8085  -0.4636  -0.1262  -0.4564
#> 
#> sigma: 0.167766
#> lambda (gene 1): 1.9112
#> epsilon (gene 1): 0.0279031
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 744.166
#> delta (gene 1): 0.605687
#> phi (cell 1): 0.876861
#> s (cell 1): 0.709601
#> nu (cell 1): 0.413092
#> theta (batch 1): 0.16913
#> betas:    0.7247  -0.1006   0.9819   0.3224   1.0603   0.2499   0.8380  -0.1893   0.1509  -0.2620   0.0183  -0.5274
#> 
#> sigma: 0.126265
#> lambda (gene 1): 0.181858
#> epsilon (gene 1): -0.24333
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 947.177
#> delta (gene 1): 0.619412
#> phi (cell 1): 0.854459
#> s (cell 1): 0.445577
#> nu (cell 1): 0.427513
#> theta (batch 1): 0.323994
#> betas:    0.6672   0.0121   0.2559   1.0516   0.7308   0.3538   0.4316  -0.2028  -0.2446  -0.1998  -0.5419  -0.4630
#> 
#> sigma: 0.197764
#> lambda (gene 1): 0.554432
#> epsilon (gene 1): -0.293731
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 592.383
#> delta (gene 1): 0.894875
#> phi (cell 1): 0.824539
#> s (cell 1): 0.236201
#> nu (cell 1): 0.423255
#> theta (batch 1): 0.286583
#> betas:    0.3126   0.0091   0.3648   1.3319   0.7966   0.5197   0.3895   0.0886   0.1840  -0.3865  -0.4537  -0.0351
#> 
#> sigma: 0.165821
#> lambda (gene 1): 1.35192
#> epsilon (gene 1): 0.103877
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 797.81
#> delta (gene 1): 0.754329
#> phi (cell 1): 0.913481
#> s (cell 1): 0.373687
#> nu (cell 1): 0.426386
#> theta (batch 1): 0.150238
#> betas:    1.3771  -0.1287   0.5023   0.3110   0.7669   0.2532   0.3264   0.0247  -0.1275  -0.0286  -0.4592  -0.1490
#> 
#> sigma: 0.106667
#> lambda (gene 1): 0.958757
#> epsilon (gene 1): -0.31487
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 700.305
#> delta (gene 1): 0.916123
#> phi (cell 1): 0.918616
#> s (cell 1): 0.342122
#> nu (cell 1): 0.427314
#> theta (batch 1): 0.117058
#> betas:    0.7735  -0.0962   0.7953   0.8284   0.4036   0.7504   0.4568   0.0850   0.0655   0.0567  -0.7743   0.4343
#> 
#> sigma: 0.1265
#> lambda (gene 1): 2.12182
#> epsilon (gene 1): 0.108017
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 632.641
#> delta (gene 1): 0.756782
#> phi (cell 1): 0.928128
#> s (cell 1): 0.46447
#> nu (cell 1): 0.422748
#> theta (batch 1): 0.00819531
#> betas:    1.0385  -0.0514   0.5558   0.8826  -0.2015   1.1345   0.0054   0.1765  -0.3538  -0.0775  -0.3634  -0.6966
#> 
#> sigma: 0.196224
#> lambda (gene 1): 1.01932
#> epsilon (gene 1): -0.228537
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 615.618
#> delta (gene 1): 0.624361
#> phi (cell 1): 0.887394
#> s (cell 1): 0.423878
#> nu (cell 1): 0.425512
#> theta (batch 1): 0.00684581
#> betas:    1.3501  -0.1398   1.1788   0.4375   0.1968   0.6581   0.1626   0.0134   0.1290  -0.6725   0.2763  -0.6893
#> 
#> sigma: 0.155231
#> lambda (gene 1): 0.503846
#> epsilon (gene 1): -0.324812
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 921.797
#> delta (gene 1): 0.928119
#> phi (cell 1): 0.754262
#> s (cell 1): 0.429749
#> nu (cell 1): 0.411987
#> theta (batch 1): 0.00258292
#> betas:    0.9753  -0.1306   1.2766   0.7011   0.3027   0.9355   0.2583   0.1369   0.0400   0.1498  -0.2906  -0.2643
#> 
#> sigma: 0.225932
#> lambda (gene 1): 0.295539
#> epsilon (gene 1): -0.0976409
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 629.403
#> delta (gene 1): 0.659262
#> phi (cell 1): 0.916312
#> s (cell 1): 0.414094
#> nu (cell 1): 0.414179
#> theta (batch 1): 0.087729
#> betas:    1.2517  -0.0650   0.2499   0.8960   0.1243   0.4726   0.1411   0.1091  -0.1813  -0.4823  -0.5232  -0.4209
#> 
#> sigma: 0.213549
#> lambda (gene 1): 1.57471
#> epsilon (gene 1): -0.145553
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 855.481
#> delta (gene 1): 0.82819
#> phi (cell 1): 0.841215
#> s (cell 1): 0.3198
#> nu (cell 1): 0.410269
#> theta (batch 1): 0.303421
#> betas:    0.9336  -0.0498   0.9964   0.7895   0.1069   0.5366   0.5556  -0.5296   0.0949  -0.2823  -0.2748  -0.6167
#> 
#> sigma: 0.157608
#> lambda (gene 1): 0.788704
#> epsilon (gene 1): -0.0908931
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 504.559
#> delta (gene 1): 0.962755
#> phi (cell 1): 0.819233
#> s (cell 1): 0.440342
#> nu (cell 1): 0.422235
#> theta (batch 1): 0.564705
#> betas:    0.5117  -0.0228   0.5248   1.0058   0.5770   0.9822   0.0240   0.2111  -0.1815  -0.0185  -0.5529  -0.0503
#> 
#> sigma: 0.154594
#> lambda (gene 1): 1.35483
#> epsilon (gene 1): 0.12647
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 889.463
#> delta (gene 1): 0.699815
#> phi (cell 1): 0.806698
#> s (cell 1): 1.02458
#> nu (cell 1): 0.42213
#> theta (batch 1): 0.430877
#> betas:    1.0076  -0.1612   0.4049   0.8475   0.2535   0.8382   0.2224   0.4250  -0.0375  -0.1637  -0.2746   0.2902
#> 
#> sigma: 0.176057
#> lambda (gene 1): 1.34146
#> epsilon (gene 1): -0.079543
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 594.147
#> delta (gene 1): 0.582771
#> phi (cell 1): 0.793596
#> s (cell 1): 0.609179
#> nu (cell 1): 0.420863
#> theta (batch 1): 0.417313
#> betas:    0.3142  -0.1083   0.9203   0.6934   0.9799   0.7209   0.9032   0.1466   0.1219  -0.0207  -0.1584   0.4606
#> 
#> sigma: 0.169638
#> lambda (gene 1): 0.233385
#> epsilon (gene 1): -0.269463
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 738.066
#> delta (gene 1): 0.973441
#> phi (cell 1): 0.792911
#> s (cell 1): 0.39347
#> nu (cell 1): 0.408744
#> theta (batch 1): 0.659388
#> betas:    0.6249  -0.0196   0.6726   0.9079   0.8293  -0.0344   0.8393  -0.2163  -0.0280  -0.1423  -0.3994  -0.4452
#> 
#> sigma: 0.152529
#> lambda (gene 1): 1.57571
#> epsilon (gene 1): 0.0896808
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 577.615
#> delta (gene 1): 0.742123
#> phi (cell 1): 0.835147
#> s (cell 1): 2.2995
#> nu (cell 1): 0.432402
#> theta (batch 1): 0.658455
#> betas:    1.0309  -0.1400   0.6770   0.3727   1.2032   0.1139   0.3460   0.2247  -0.0298  -0.1059  -0.3148  -0.0751
#> 
#> sigma: 0.144354
#> lambda (gene 1): 1.45493
#> epsilon (gene 1): -0.070939
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 725.978
#> delta (gene 1): 0.943316
#> phi (cell 1): 0.818702
#> s (cell 1): 0.42897
#> nu (cell 1): 0.432789
#> theta (batch 1): 0.668368
#> betas:    0.6513  -0.0704   0.2069   0.8370   0.6459   0.3406   0.9143  -0.2372   0.0764   0.2770  -0.6363  -0.3040
#> 
#> sigma: 0.165549
#> lambda (gene 1): 3.00449
#> epsilon (gene 1): 0.0251229
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 714.003
#> delta (gene 1): 0.893468
#> phi (cell 1): 0.697605
#> s (cell 1): 0.312047
#> nu (cell 1): 0.42799
#> theta (batch 1): 0.419937
#> betas:    1.2443  -0.1637   0.7849   0.6943   0.9438   0.0577   0.7196  -0.0574   0.1138  -0.2810  -0.2836   0.0133
#> 
#> sigma: 0.144329
#> lambda (gene 1): 1.2546
#> epsilon (gene 1): 0.138027
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 497.35
#> delta (gene 1): 0.640009
#> phi (cell 1): 0.709636
#> s (cell 1): 0.194999
#> nu (cell 1): 0.422915
#> theta (batch 1): 0.350329
#> betas:    1.1578  -0.1767   0.0681   1.1435   0.2547   0.7241   0.0861   0.1748   0.5950  -0.7076  -0.0594   0.1241
#> 
#> sigma: 0.165198
#> lambda (gene 1): 1.10882
#> epsilon (gene 1): -0.142664
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 586.075
#> delta (gene 1): 0.715551
#> phi (cell 1): 0.765597
#> s (cell 1): 0.581674
#> nu (cell 1): 0.413508
#> theta (batch 1): 0.293895
#> betas:    0.4742  -0.0396   1.2816   1.1740   0.3903   0.5845   0.2804   0.2358   0.3158  -0.7720   0.0676  -0.2268
#> 
#> sigma: 0.170395
#> lambda (gene 1): 0.407918
#> epsilon (gene 1): -0.0285883
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 753.955
#> delta (gene 1): 0.738074
#> phi (cell 1): 0.741529
#> s (cell 1): 0.397455
#> nu (cell 1): 0.416191
#> theta (batch 1): 0.408387
#> betas:    0.5754   0.0005   0.3671   1.5614   0.5014   0.4416   0.1301   0.2907  -0.1212  -0.5144   0.1647  -0.6216
#> 
#> sigma: 0.270966
#> lambda (gene 1): 1.17974
#> epsilon (gene 1): -0.303176
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 14980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 705.473
#> delta (gene 1): 0.65927
#> phi (cell 1): 0.787821
#> s (cell 1): 0.715489
#> nu (cell 1): 0.413041
#> theta (batch 1): 0.334791
#> betas:    0.4388  -0.0699   0.8858   1.1196   0.3115   0.7505   0.5663   0.1224   0.2321  -0.4008   0.0558  -0.1744
#> 
#> sigma: 0.152495
#> lambda (gene 1): 1.70543
#> epsilon (gene 1): -0.181103
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 555.31
#> delta (gene 1): 0.678109
#> phi (cell 1): 0.844458
#> s (cell 1): 0.632795
#> nu (cell 1): 0.424068
#> theta (batch 1): 0.2736
#> betas:    0.5757   0.0382   0.7866   0.7602   0.9614   0.4490  -0.0234   0.1680  -0.1178  -0.4558  -0.4548  -0.5312
#> 
#> sigma: 0.178018
#> lambda (gene 1): 0.754529
#> epsilon (gene 1): -0.22266
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 719.487
#> delta (gene 1): 0.771337
#> phi (cell 1): 0.798033
#> s (cell 1): 0.332166
#> nu (cell 1): 0.419443
#> theta (batch 1): 0.448921
#> betas:    0.8416  -0.1025   0.9731   0.4650   0.5599   0.5859   0.7936  -0.2545   0.0964  -0.2618  -0.3053   0.2197
#> 
#> sigma: 0.134533
#> lambda (gene 1): 0.320329
#> epsilon (gene 1): -0.0274032
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 555.996
#> delta (gene 1): 1.02821
#> phi (cell 1): 0.837344
#> s (cell 1): 0.283773
#> nu (cell 1): 0.419748
#> theta (batch 1): 0.297311
#> betas:    0.1833   0.1036   0.7177   0.8960   0.8655   0.5330   0.4723  -0.1491  -0.2280  -0.2855  -0.7356  -0.5574
#> 
#> sigma: 0.202035
#> lambda (gene 1): 0.91132
#> epsilon (gene 1): 0.324609
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 713.897
#> delta (gene 1): 0.818418
#> phi (cell 1): 0.786456
#> s (cell 1): 1.4822
#> nu (cell 1): 0.424285
#> theta (batch 1): 0.264715
#> betas:    0.5913   0.0458   0.3034   0.9640   0.7090   0.2713  -0.0126   0.1380  -0.1398  -0.3398  -0.8258  -0.6638
#> 
#> sigma: 0.16687
#> lambda (gene 1): 2.2439
#> epsilon (gene 1): 0.126519
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 1026.59
#> delta (gene 1): 0.59179
#> phi (cell 1): 0.801528
#> s (cell 1): 0.709813
#> nu (cell 1): 0.417031
#> theta (batch 1): 0.15107
#> betas:    0.7240  -0.0415   0.7972   0.3975   0.6716   0.6103   0.1738   0.6675  -0.6327  -0.0261  -0.4739  -0.1571
#> 
#> sigma: 0.175583
#> lambda (gene 1): 2.08418
#> epsilon (gene 1): -0.277112
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 562.1
#> delta (gene 1): 0.872174
#> phi (cell 1): 0.803547
#> s (cell 1): 2.69758
#> nu (cell 1): 0.425455
#> theta (batch 1): 0.445519
#> betas:    1.0623  -0.0521   0.7301   0.7561   0.3174   0.8511  -0.4051   0.5760  -0.6087  -0.0547  -0.5568  -0.5706
#> 
#> sigma: 0.130327
#> lambda (gene 1): 3.24157
#> epsilon (gene 1): 0.0781875
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 689.391
#> delta (gene 1): 0.751753
#> phi (cell 1): 0.820584
#> s (cell 1): 0.513854
#> nu (cell 1): 0.416334
#> theta (batch 1): 0.0763027
#> betas:    0.5618  -0.0261   0.0035   1.3121   0.4243   0.9291  -0.0163   0.2576   0.0106  -0.3914  -0.3521  -0.3178
#> 
#> sigma: 0.200721
#> lambda (gene 1): 2.21875
#> epsilon (gene 1): 0.0140393
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 742.079
#> delta (gene 1): 0.910525
#> phi (cell 1): 0.803777
#> s (cell 1): 0.451196
#> nu (cell 1): 0.417733
#> theta (batch 1): 0.213911
#> betas:    1.0281  -0.1184   0.4992   1.2076   0.2352   0.8314   0.1034   0.0212   0.2391  -0.1662  -0.3152  -0.4291
#> 
#> sigma: 0.169366
#> lambda (gene 1): 0.865871
#> epsilon (gene 1): 0.0313577
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 648.123
#> delta (gene 1): 0.726995
#> phi (cell 1): 0.783798
#> s (cell 1): 0.200374
#> nu (cell 1): 0.427209
#> theta (batch 1): 0.308729
#> betas:    2.0937  -0.0585   0.5125   0.4566  -0.1655   0.4102  -0.4187  -0.5001  -0.2155  -0.6886  -0.9166  -0.5557
#> 
#> sigma: 0.237032
#> lambda (gene 1): 0.167034
#> epsilon (gene 1): -0.218108
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 610.911
#> delta (gene 1): 0.653261
#> phi (cell 1): 0.841581
#> s (cell 1): 0.320155
#> nu (cell 1): 0.423524
#> theta (batch 1): 0.532325
#> betas:    0.5984   0.0552   1.2852   0.1138   1.0975   0.3663  -0.0767   0.0499  -0.0841  -0.3076  -0.5626  -0.8471
#> 
#> sigma: 0.116804
#> lambda (gene 1): 1.3
#> epsilon (gene 1): -0.363365
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 666.292
#> delta (gene 1): 0.718768
#> phi (cell 1): 0.856399
#> s (cell 1): 0.266642
#> nu (cell 1): 0.405113
#> theta (batch 1): 0.44523
#> betas:    0.5479  -0.0533   1.2055   0.8770   0.6865   0.5333   0.5710   0.1783   0.0575  -0.3698  -0.0762  -0.4448
#> 
#> sigma: 0.141967
#> lambda (gene 1): 1.11124
#> epsilon (gene 1): -0.064053
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 564.573
#> delta (gene 1): 0.464986
#> phi (cell 1): 0.755963
#> s (cell 1): 0.744767
#> nu (cell 1): 0.443796
#> theta (batch 1): 0.175271
#> betas:    1.1367  -0.1456   0.4602   0.8021   0.6052   0.7805   0.1774   0.0394   0.0817  -0.1647  -0.4941   0.0883
#> 
#> sigma: 0.156387
#> lambda (gene 1): 0.308405
#> epsilon (gene 1): -0.488279
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 559.584
#> delta (gene 1): 0.691387
#> phi (cell 1): 0.805005
#> s (cell 1): 0.455848
#> nu (cell 1): 0.424559
#> theta (batch 1): 0.0505925
#> betas:    1.6738  -0.1999   0.4635   0.8110   0.4713   0.1770   0.3463   0.0099   0.0773  -0.6161  -0.0622  -0.3860
#> 
#> sigma: 0.119619
#> lambda (gene 1): 0.998961
#> epsilon (gene 1): -0.00422949
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 636.027
#> delta (gene 1): 0.692643
#> phi (cell 1): 0.685428
#> s (cell 1): 0.529
#> nu (cell 1): 0.435962
#> theta (batch 1): 0.0386993
#> betas:    1.2205  -0.0614   1.0644   0.5840   0.2455  -0.0463   0.8754  -0.5921   0.1032  -0.5800  -0.2941  -0.6434
#> 
#> sigma: 0.167053
#> lambda (gene 1): 0.875668
#> epsilon (gene 1): -0.162005
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 630.862
#> delta (gene 1): 0.684959
#> phi (cell 1): 0.673791
#> s (cell 1): 0.519596
#> nu (cell 1): 0.423544
#> theta (batch 1): 0.0189579
#> betas:    0.1745   0.0936   0.5068   0.8710   0.9905   0.5397   0.1452   0.0517   0.0344  -0.2191  -0.8080  -0.5978
#> 
#> sigma: 0.125006
#> lambda (gene 1): 1.18407
#> epsilon (gene 1): -0.207003
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 687.707
#> delta (gene 1): 0.91156
#> phi (cell 1): 0.710617
#> s (cell 1): 0.528644
#> nu (cell 1): 0.438215
#> theta (batch 1): 0.0116283
#> betas:    0.1961   0.0414   0.4070   1.0307   1.2771   0.1310   0.1248   0.4321  -0.1904  -0.3821  -0.0444  -0.7812
#> 
#> sigma: 0.184618
#> lambda (gene 1): 1.87482
#> epsilon (gene 1): 0.105561
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 627.397
#> delta (gene 1): 0.833226
#> phi (cell 1): 0.868293
#> s (cell 1): 0.493762
#> nu (cell 1): 0.41429
#> theta (batch 1): 0.0217488
#> betas:    1.1920  -0.1165   0.0280   1.1662   0.6692   0.2455   0.2284   0.1183  -0.0369  -0.2947  -0.4393  -0.1243
#> 
#> sigma: 0.132322
#> lambda (gene 1): 1.56692
#> epsilon (gene 1): 0.0558107
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 733.49
#> delta (gene 1): 1.0496
#> phi (cell 1): 0.793283
#> s (cell 1): 0.420766
#> nu (cell 1): 0.414857
#> theta (batch 1): 0.011726
#> betas:    0.8234  -0.0843   0.0406   0.9960   0.5628   0.5265   0.4727  -0.0743   0.2024  -0.2468  -0.4425  -0.1003
#> 
#> sigma: 0.183691
#> lambda (gene 1): 0.359746
#> epsilon (gene 1): 0.264883
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 619.337
#> delta (gene 1): 0.757347
#> phi (cell 1): 0.687813
#> s (cell 1): 0.374504
#> nu (cell 1): 0.428542
#> theta (batch 1): 0.0119357
#> betas:    0.9520   0.1358   0.6549   0.8678   0.1239   0.5388   0.1540  -1.1899   0.5451  -0.6317  -1.3172  -0.8926
#> 
#> sigma: 0.218078
#> lambda (gene 1): 0.499119
#> epsilon (gene 1): -0.348827
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 683.32
#> delta (gene 1): 0.87307
#> phi (cell 1): 0.67927
#> s (cell 1): 0.448715
#> nu (cell 1): 0.441714
#> theta (batch 1): 0.00808839
#> betas:    1.1033  -0.1035   0.9079  -0.0076   1.1297  -0.0330   0.7625  -0.0814  -0.1240  -0.3909  -0.0189  -0.7655
#> 
#> sigma: 0.165939
#> lambda (gene 1): 0.729187
#> epsilon (gene 1): 0.147641
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 654.264
#> delta (gene 1): 0.903346
#> phi (cell 1): 0.77494
#> s (cell 1): 0.373347
#> nu (cell 1): 0.420132
#> theta (batch 1): 0.0299493
#> betas:    0.2490  -0.0067   1.0894   0.8329   0.6341   0.8684   0.2119   0.1010   0.2935  -0.3950   0.0717  -0.6027
#> 
#> sigma: 0.138011
#> lambda (gene 1): 0.869483
#> epsilon (gene 1): -0.0132514
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 689.636
#> delta (gene 1): 0.616166
#> phi (cell 1): 0.777923
#> s (cell 1): 0.393554
#> nu (cell 1): 0.408598
#> theta (batch 1): 0.0533282
#> betas:    0.8452  -0.0744   0.1911   1.0607   0.4904   0.2162   0.6045  -0.0426   0.0723  -0.4655  -0.4726   0.1486
#> 
#> sigma: 0.155831
#> lambda (gene 1): 1.86932
#> epsilon (gene 1): -0.113615
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 742.927
#> delta (gene 1): 0.744674
#> phi (cell 1): 0.72771
#> s (cell 1): 0.269878
#> nu (cell 1): 0.424811
#> theta (batch 1): 0.246519
#> betas:    0.9299  -0.0059   0.2248   1.0045   0.2756   0.5585   0.2083  -0.1656  -0.1470  -0.5796  -0.4568  -0.3651
#> 
#> sigma: 0.14463
#> lambda (gene 1): 2.0254
#> epsilon (gene 1): -0.0319382
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 739.158
#> delta (gene 1): 0.556252
#> phi (cell 1): 0.749464
#> s (cell 1): 0.766822
#> nu (cell 1): 0.430815
#> theta (batch 1): 0.332013
#> betas:    1.4212  -0.2056   0.2099   0.3130   0.7775   0.7833  -0.0501   0.2991  -0.0042  -0.0856  -0.0632  -0.1210
#> 
#> sigma: 0.17614
#> lambda (gene 1): 0.660298
#> epsilon (gene 1): -0.532313
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 690.922
#> delta (gene 1): 0.623474
#> phi (cell 1): 0.765048
#> s (cell 1): 0.56681
#> nu (cell 1): 0.412928
#> theta (batch 1): 0.310807
#> betas:    1.3095  -0.2884   0.8476   0.2114   0.8406   0.7288   0.2394   0.4205   0.2475   0.0478   0.1644   0.0986
#> 
#> sigma: 0.15491
#> lambda (gene 1): 1.26564
#> epsilon (gene 1): -0.324586
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 639.166
#> delta (gene 1): 0.940448
#> phi (cell 1): 0.783198
#> s (cell 1): 0.486564
#> nu (cell 1): 0.408021
#> theta (batch 1): 0.147406
#> betas:    0.8664  -0.1111   0.7079   0.7019   0.6161   0.6866   0.0859   0.4937  -0.4359   0.0400  -0.3126  -0.3107
#> 
#> sigma: 0.136231
#> lambda (gene 1): 0.837235
#> epsilon (gene 1): 0.274593
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 558.41
#> delta (gene 1): 0.670748
#> phi (cell 1): 0.863684
#> s (cell 1): 0.874819
#> nu (cell 1): 0.425045
#> theta (batch 1): 0.34779
#> betas:    0.3278   0.0482   0.6175   0.4299   1.3445   0.2842   0.3584  -0.0840   0.1234  -0.0322  -1.0909  -0.1568
#> 
#> sigma: 0.175962
#> lambda (gene 1): 0.0986349
#> epsilon (gene 1): -0.213569
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 670.264
#> delta (gene 1): 1.08229
#> phi (cell 1): 0.765809
#> s (cell 1): 0.375202
#> nu (cell 1): 0.432051
#> theta (batch 1): 0.318819
#> betas:    1.1445  -0.2361  -0.1196   0.9482   0.3556   1.1364  -0.0374   0.6350   0.1497  -0.1045   0.0217   0.2116
#> 
#> sigma: 0.186055
#> lambda (gene 1): 1.80609
#> epsilon (gene 1): 0.298394
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 811.014
#> delta (gene 1): 0.573563
#> phi (cell 1): 0.736434
#> s (cell 1): 0.709578
#> nu (cell 1): 0.428018
#> theta (batch 1): 0.0409461
#> betas:    0.9209  -0.0623   0.7681   0.4557   0.2009   1.0916   0.1032   0.1128  -0.2190  -0.2763   0.0135  -0.8235
#> 
#> sigma: 0.169432
#> lambda (gene 1): 0.780909
#> epsilon (gene 1): -0.44583
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 650.729
#> delta (gene 1): 0.890201
#> phi (cell 1): 0.66204
#> s (cell 1): 0.40724
#> nu (cell 1): 0.437704
#> theta (batch 1): 0.0452682
#> betas:    0.4135  -0.0145   1.1173   0.7807   0.8237   0.2469   0.7542  -0.0366  -0.3645   0.5183  -0.8085  -0.2265
#> 
#> sigma: 0.196501
#> lambda (gene 1): 1.42548
#> epsilon (gene 1): -0.0492461
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 817.766
#> delta (gene 1): 0.642659
#> phi (cell 1): 0.817401
#> s (cell 1): 0.391794
#> nu (cell 1): 0.426214
#> theta (batch 1): 0.0436883
#> betas:    1.1339  -0.0593   0.5096   0.7546   0.4544   0.4966   0.0760   0.1526  -0.3524  -0.3510  -0.5879  -0.2434
#> 
#> sigma: 0.12726
#> lambda (gene 1): 2.02469
#> epsilon (gene 1): -0.120363
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 696.051
#> delta (gene 1): 0.712786
#> phi (cell 1): 0.83826
#> s (cell 1): 0.440338
#> nu (cell 1): 0.405209
#> theta (batch 1): 0.0163429
#> betas:    0.5698  -0.0108   1.0006   0.9923   0.8751   0.4963   0.1016   0.5663  -0.4997   0.0229  -0.2144  -0.9681
#> 
#> sigma: 0.172243
#> lambda (gene 1): 0.558912
#> epsilon (gene 1): -0.213937
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 763.664
#> delta (gene 1): 0.696775
#> phi (cell 1): 0.907001
#> s (cell 1): 0.440104
#> nu (cell 1): 0.429214
#> theta (batch 1): 0.00456894
#> betas:    0.9192  -0.1919   0.5770   0.8874   0.5139   0.7720   0.3696   0.6759  -0.2455   0.0975  -0.0425   0.1098
#> 
#> sigma: 0.131416
#> lambda (gene 1): 0.695943
#> epsilon (gene 1): -0.0986857
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 680.573
#> delta (gene 1): 0.599335
#> phi (cell 1): 0.92232
#> s (cell 1): 0.467965
#> nu (cell 1): 0.421835
#> theta (batch 1): 0.0147735
#> betas:    1.0422  -0.2594   1.1852   0.5889   0.6916   0.4112   1.1066  -0.1984   0.6576  -0.1638   0.3429   0.1970
#> 
#> sigma: 0.181699
#> lambda (gene 1): 1.54945
#> epsilon (gene 1): -0.407989
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 666.816
#> delta (gene 1): 0.826257
#> phi (cell 1): 0.861388
#> s (cell 1): 0.418902
#> nu (cell 1): 0.414069
#> theta (batch 1): 0.0319227
#> betas:    0.5479  -0.1534   0.8366   0.7728   0.6775   0.8577   0.7329   0.5483  -0.3873   0.4023   0.1841  -0.4103
#> 
#> sigma: 0.198359
#> lambda (gene 1): 1.34119
#> epsilon (gene 1): -0.0309974
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 820.675
#> delta (gene 1): 0.647314
#> phi (cell 1): 0.770763
#> s (cell 1): 0.506113
#> nu (cell 1): 0.423498
#> theta (batch 1): 0.0455412
#> betas:    0.8625  -0.0560   0.7964   0.8689   0.5780   0.5482   0.1883   0.4162  -0.3684  -0.2277  -0.3457  -0.5728
#> 
#> sigma: 0.132114
#> lambda (gene 1): 1.08554
#> epsilon (gene 1): -0.124233
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 696.328
#> delta (gene 1): 0.530741
#> phi (cell 1): 0.831244
#> s (cell 1): 0.48869
#> nu (cell 1): 0.435364
#> theta (batch 1): 0.345468
#> betas:    1.0040  -0.0150   0.0758   0.2924   0.5041   0.7146   0.0501   0.2137  -0.4710  -0.3933  -0.4621  -0.6949
#> 
#> sigma: 0.139208
#> lambda (gene 1): 0.375661
#> epsilon (gene 1): -0.334181
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 710.787
#> delta (gene 1): 0.643936
#> phi (cell 1): 0.835858
#> s (cell 1): 1.14473
#> nu (cell 1): 0.415574
#> theta (batch 1): 0.373669
#> betas:    0.9841   0.0255   1.1254   0.2910   0.5411   0.2842  -0.0430   0.1164  -0.2175  -0.8155  -0.4090  -0.6912
#> 
#> sigma: 0.146286
#> lambda (gene 1): 0.919342
#> epsilon (gene 1): -0.151338
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 752.406
#> delta (gene 1): 0.973897
#> phi (cell 1): 0.779555
#> s (cell 1): 0.426447
#> nu (cell 1): 0.430394
#> theta (batch 1): 0.154287
#> betas:    0.6716  -0.0936   0.4584   1.0143   0.6143   0.4487   0.3968   0.5318  -0.0477  -0.4025   0.3514  -0.5447
#> 
#> sigma: 0.180211
#> lambda (gene 1): 1.01974
#> epsilon (gene 1): 0.119882
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 904.975
#> delta (gene 1): 1.06166
#> phi (cell 1): 0.738941
#> s (cell 1): 0.5151
#> nu (cell 1): 0.425562
#> theta (batch 1): 0.0959608
#> betas:    1.2280  -0.1469   0.3143   1.1316   0.4591  -0.0476   0.7285   0.2460  -0.3980  -0.2828  -0.0091  -0.3239
#> 
#> sigma: 0.163562
#> lambda (gene 1): 0.778761
#> epsilon (gene 1): 0.373623
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 713.858
#> delta (gene 1): 0.730706
#> phi (cell 1): 0.800513
#> s (cell 1): 0.312975
#> nu (cell 1): 0.418593
#> theta (batch 1): 0.0937625
#> betas:    0.8909  -0.0774   1.1677   0.6450   0.4127   0.6426   0.5345  -0.1429  -0.0915   0.0882  -0.6070  -0.3332
#> 
#> sigma: 0.180608
#> lambda (gene 1): 0.875363
#> epsilon (gene 1): -0.148048
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 575.69
#> delta (gene 1): 0.90383
#> phi (cell 1): 0.853281
#> s (cell 1): 0.400863
#> nu (cell 1): 0.412413
#> theta (batch 1): 0.263566
#> betas:    1.2385  -0.2248   0.0697   0.6322   0.5807   0.7406   0.2469   0.3716  -0.4866   0.6327  -0.1890  -0.3211
#> 
#> sigma: 0.174926
#> lambda (gene 1): 0.375628
#> epsilon (gene 1): -0.00927665
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 641.782
#> delta (gene 1): 0.750737
#> phi (cell 1): 0.710604
#> s (cell 1): 0.446025
#> nu (cell 1): 0.426908
#> theta (batch 1): 0.374801
#> betas:    1.3010  -0.1055   0.6611  -0.0181   1.0540   0.1864  -0.1149   0.1598   0.2119  -0.4673  -0.5699  -0.3928
#> 
#> sigma: 0.160586
#> lambda (gene 1): 0.931936
#> epsilon (gene 1): -0.0394556
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 833.715
#> delta (gene 1): 0.73033
#> phi (cell 1): 0.691826
#> s (cell 1): 0.378263
#> nu (cell 1): 0.426304
#> theta (batch 1): 0.198956
#> betas:    1.0146  -0.0989   0.9063   0.2498   0.8730   0.3725   0.4683  -0.4668   0.3575  -0.5716  -0.1148  -0.2115
#> 
#> sigma: 0.200446
#> lambda (gene 1): 1.05674
#> epsilon (gene 1): -0.103148
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 574.531
#> delta (gene 1): 1.1702
#> phi (cell 1): 0.779548
#> s (cell 1): 0.639041
#> nu (cell 1): 0.436353
#> theta (batch 1): 0.272407
#> betas:    0.6348  -0.0829   0.1303   0.9278   0.6658   0.7468   0.4135  -0.0808   0.3282  -0.2468   0.0150  -0.4015
#> 
#> sigma: 0.163407
#> lambda (gene 1): 0.83308
#> epsilon (gene 1): 0.212774
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 618.881
#> delta (gene 1): 0.632819
#> phi (cell 1): 0.802446
#> s (cell 1): 0.711017
#> nu (cell 1): 0.414276
#> theta (batch 1): 0.158207
#> betas:    1.2021  -0.0741   0.8897   0.3377   0.5093   0.8243  -0.1784  -0.0256   0.0461  -0.4245  -0.5603  -0.2574
#> 
#> sigma: 0.161038
#> lambda (gene 1): 1.31737
#> epsilon (gene 1): -0.262624
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 885.877
#> delta (gene 1): 0.547628
#> phi (cell 1): 0.679452
#> s (cell 1): 0.612718
#> nu (cell 1): 0.438742
#> theta (batch 1): 0.126354
#> betas:    0.8354  -0.0838   1.0708   0.6883   0.8322   0.0950   0.5754  -0.3061   0.4566  -0.0913  -0.5274  -0.2201
#> 
#> sigma: 0.195624
#> lambda (gene 1): 0.446451
#> epsilon (gene 1): -0.574437
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 689.519
#> delta (gene 1): 0.695382
#> phi (cell 1): 0.766546
#> s (cell 1): 0.325601
#> nu (cell 1): 0.432346
#> theta (batch 1): 0.189159
#> betas:    0.7816  -0.0301   0.4220   0.7807   0.4598   0.6503   0.1108   0.0934   0.0016  -0.5510  -0.2678  -0.4835
#> 
#> sigma: 0.108557
#> lambda (gene 1): 2.98972
#> epsilon (gene 1): -0.0769514
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 750.047
#> delta (gene 1): 0.750594
#> phi (cell 1): 0.75115
#> s (cell 1): 0.312915
#> nu (cell 1): 0.415022
#> theta (batch 1): 0.292385
#> betas:    0.5288  -0.1328   1.3368   0.9071   0.8335   0.6054   0.5188   0.4110   0.2228  -0.1419  -0.0262   0.1302
#> 
#> sigma: 0.184712
#> lambda (gene 1): 2.7082
#> epsilon (gene 1): -0.0593757
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 15980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 573.4
#> delta (gene 1): 0.863513
#> phi (cell 1): 0.799532
#> s (cell 1): 0.438015
#> nu (cell 1): 0.428986
#> theta (batch 1): 0.254066
#> betas:    0.9333  -0.0774   0.7627   0.1948   0.9980   0.3019   0.3037  -0.4511   0.7778  -0.6949  -0.4683  -0.0662
#> 
#> sigma: 0.127074
#> lambda (gene 1): 3.2957
#> epsilon (gene 1): 0.0825007
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 785.034
#> delta (gene 1): 0.776245
#> phi (cell 1): 0.776147
#> s (cell 1): 0.246975
#> nu (cell 1): 0.421829
#> theta (batch 1): 0.148097
#> betas:    0.9563  -0.1126   0.1000   0.8062   0.6012   0.6793   0.3622  -0.0746   0.2454  -0.3739  -0.1712  -0.0667
#> 
#> sigma: 0.187489
#> lambda (gene 1): 0.928626
#> epsilon (gene 1): -0.0959246
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 617.185
#> delta (gene 1): 1.14285
#> phi (cell 1): 0.770862
#> s (cell 1): 0.17198
#> nu (cell 1): 0.427162
#> theta (batch 1): 0.267318
#> betas:    0.5916  -0.0777   0.6948   0.6956   0.8967   0.6232   0.4233  -0.1843   0.3928  -0.3135  -0.2413   0.0854
#> 
#> sigma: 0.116209
#> lambda (gene 1): 0.721556
#> epsilon (gene 1): 0.297498
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 667.796
#> delta (gene 1): 0.587363
#> phi (cell 1): 0.822265
#> s (cell 1): 0.479644
#> nu (cell 1): 0.411142
#> theta (batch 1): 0.0626132
#> betas:    0.8641  -0.0517   0.6665   0.5961   0.6126   0.7089   0.2459   0.2248  -0.2218  -0.2135  -0.5387  -0.3390
#> 
#> sigma: 0.126395
#> lambda (gene 1): 1.03211
#> epsilon (gene 1): -0.244266
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 801.854
#> delta (gene 1): 0.636045
#> phi (cell 1): 0.882065
#> s (cell 1): 0.546111
#> nu (cell 1): 0.422566
#> theta (batch 1): 0.0518485
#> betas:    1.3187   0.0211   0.4602   0.4774   0.3849   0.2618   0.4262  -0.5437  -0.3407  -0.3536  -1.1365  -0.2518
#> 
#> sigma: 0.165906
#> lambda (gene 1): 0.879699
#> epsilon (gene 1): -0.299261
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 667.062
#> delta (gene 1): 0.726323
#> phi (cell 1): 0.950391
#> s (cell 1): 0.495764
#> nu (cell 1): 0.414653
#> theta (batch 1): 0.07741
#> betas:    0.6541   0.0474   0.7035   0.6531   0.7705   0.6080  -0.1852   0.1118  -0.2585  -0.0096  -0.8713  -0.9113
#> 
#> sigma: 0.163248
#> lambda (gene 1): 1.26429
#> epsilon (gene 1): -0.200985
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 682.434
#> delta (gene 1): 0.820916
#> phi (cell 1): 0.939085
#> s (cell 1): 0.729561
#> nu (cell 1): 0.427362
#> theta (batch 1): 0.110318
#> betas:    1.1525  -0.1157   0.0883   0.6142   0.6495   0.3742   0.4820  -0.1348  -0.0848  -0.2419  -0.2834  -0.4100
#> 
#> sigma: 0.109029
#> lambda (gene 1): 0.815369
#> epsilon (gene 1): 0.0626961
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 699.889
#> delta (gene 1): 0.696868
#> phi (cell 1): 0.842416
#> s (cell 1): 0.502814
#> nu (cell 1): 0.411347
#> theta (batch 1): 0.152056
#> betas:    1.2423  -0.0750   0.9830   0.2314   0.6010   0.1876   0.3563  -0.5964   0.3833  -0.4836  -0.2933  -0.9200
#> 
#> sigma: 0.144189
#> lambda (gene 1): 0.920999
#> epsilon (gene 1): -0.266807
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 609.763
#> delta (gene 1): 0.605304
#> phi (cell 1): 0.905192
#> s (cell 1): 0.513409
#> nu (cell 1): 0.408687
#> theta (batch 1): 0.0408714
#> betas:    1.0414  -0.1334   0.9086   0.8687   0.2729   0.5411   0.4380   0.2125   0.0396  -0.2697  -0.3270   0.0297
#> 
#> sigma: 0.126093
#> lambda (gene 1): 0.503372
#> epsilon (gene 1): -0.233487
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 736.993
#> delta (gene 1): 0.656005
#> phi (cell 1): 0.792954
#> s (cell 1): 0.335445
#> nu (cell 1): 0.409348
#> theta (batch 1): 0.0361227
#> betas:    0.4800  -0.1483   0.6871   0.5706   1.4299   0.3744   0.4607   0.5276   0.0022   0.2549  -0.0841   0.1032
#> 
#> sigma: 0.191199
#> lambda (gene 1): 1.93969
#> epsilon (gene 1): -0.264242
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 582.256
#> delta (gene 1): 0.60836
#> phi (cell 1): 0.698421
#> s (cell 1): 0.391055
#> nu (cell 1): 0.42364
#> theta (batch 1): 0.101099
#> betas:    0.5755  -0.0973   0.4946   0.6756   0.6695   0.4734   0.6789   0.4197  -0.3975   0.3249  -0.3661   0.2234
#> 
#> sigma: 0.138222
#> lambda (gene 1): 1.17527
#> epsilon (gene 1): -0.390155
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 678.343
#> delta (gene 1): 0.545848
#> phi (cell 1): 0.670761
#> s (cell 1): 0.178897
#> nu (cell 1): 0.434506
#> theta (batch 1): 0.248358
#> betas:    1.1355  -0.1935   0.4082   0.4544   0.9300   0.4734   0.3834   0.2836  -0.0812  -0.0277  -0.2234   0.2046
#> 
#> sigma: 0.176138
#> lambda (gene 1): 0.754611
#> epsilon (gene 1): -0.337103
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 697.694
#> delta (gene 1): 0.581588
#> phi (cell 1): 0.669181
#> s (cell 1): 0.558851
#> nu (cell 1): 0.437278
#> theta (batch 1): 0.270152
#> betas:    0.2021  -0.0709   0.6510   1.4831   0.6932   0.6413   0.9649  -0.3651   0.7844  -0.4570   0.0971  -0.1851
#> 
#> sigma: 0.16953
#> lambda (gene 1): 0.457956
#> epsilon (gene 1): -0.3042
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 550.816
#> delta (gene 1): 0.685659
#> phi (cell 1): 0.644558
#> s (cell 1): 0.395518
#> nu (cell 1): 0.427717
#> theta (batch 1): 0.0337943
#> betas:    0.5390  -0.0221  -0.0667   1.2069   0.5624   1.1157  -0.2475  -0.0042   0.7352  -1.1175   0.0262  -0.2104
#> 
#> sigma: 0.155413
#> lambda (gene 1): 1.16177
#> epsilon (gene 1): -0.0866678
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 726.262
#> delta (gene 1): 0.625319
#> phi (cell 1): 0.643227
#> s (cell 1): 0.400716
#> nu (cell 1): 0.420724
#> theta (batch 1): 0.0433578
#> betas:    0.9137   0.0038   0.5542   0.6517   0.5468   0.8169  -0.6452   0.4122  -0.0973  -0.4060  -0.8434  -0.5061
#> 
#> sigma: 0.13057
#> lambda (gene 1): 2.14035
#> epsilon (gene 1): -0.217424
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 661.975
#> delta (gene 1): 0.566654
#> phi (cell 1): 0.7217
#> s (cell 1): 0.469021
#> nu (cell 1): 0.423168
#> theta (batch 1): 0.056913
#> betas:    1.0200  -0.1060   0.3330   1.0455   0.5358   0.2603   0.2166   0.0687   0.2864  -0.4742  -0.3574   0.0401
#> 
#> sigma: 0.132428
#> lambda (gene 1): 1.17381
#> epsilon (gene 1): -0.358208
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 658.388
#> delta (gene 1): 0.567339
#> phi (cell 1): 0.724674
#> s (cell 1): 0.489045
#> nu (cell 1): 0.433919
#> theta (batch 1): 0.011267
#> betas:    0.6869  -0.0068   0.9343   0.2119   0.7395   0.8389  -0.1941   0.1428  -0.1928  -0.1971  -0.6346  -0.5359
#> 
#> sigma: 0.140564
#> lambda (gene 1): 2.52776
#> epsilon (gene 1): -0.27296
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 601.457
#> delta (gene 1): 0.828709
#> phi (cell 1): 0.817529
#> s (cell 1): 0.415404
#> nu (cell 1): 0.416574
#> theta (batch 1): 0.0582005
#> betas:    1.4285  -0.1796   0.7394   0.7572  -0.1585   0.8863   0.5145  -0.3563   0.3364  -0.5573   0.0573  -0.2194
#> 
#> sigma: 0.170507
#> lambda (gene 1): 0.321901
#> epsilon (gene 1): -0.00928893
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 770.344
#> delta (gene 1): 0.496774
#> phi (cell 1): 0.877467
#> s (cell 1): 0.429616
#> nu (cell 1): 0.41325
#> theta (batch 1): 0.0969568
#> betas:    0.8419  -0.0690   0.5802   0.8954   1.0341  -0.0303   0.9164  -0.5312   0.1738  -0.2730  -0.4822  -0.0799
#> 
#> sigma: 0.18511
#> lambda (gene 1): 1.32751
#> epsilon (gene 1): -0.451549
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 569.759
#> delta (gene 1): 1.07631
#> phi (cell 1): 0.834764
#> s (cell 1): 0.374613
#> nu (cell 1): 0.419351
#> theta (batch 1): 0.204877
#> betas:    0.8579  -0.0679   0.4695   0.6698   0.3647   0.5408   0.6527  -0.3540   0.2968  -0.5396  -0.0412  -0.6670
#> 
#> sigma: 0.156706
#> lambda (gene 1): 1.01691
#> epsilon (gene 1): 0.288992
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 698.5
#> delta (gene 1): 1.08867
#> phi (cell 1): 0.836616
#> s (cell 1): 0.463105
#> nu (cell 1): 0.422613
#> theta (batch 1): 0.266412
#> betas:    0.5197  -0.0665   0.5297   0.8927   0.9602   0.7447   0.2192   0.4187  -0.1386  -0.1919  -0.3241   0.0436
#> 
#> sigma: 0.139212
#> lambda (gene 1): 0.9848
#> epsilon (gene 1): 0.424651
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 684.803
#> delta (gene 1): 0.530284
#> phi (cell 1): 0.892793
#> s (cell 1): 1.71913
#> nu (cell 1): 0.430789
#> theta (batch 1): 0.18835
#> betas:    0.3333   0.1250   0.8044   0.9245   0.1680   0.7851   0.3851  -0.3785  -0.1798  -0.3799  -0.9845  -0.4347
#> 
#> sigma: 0.102995
#> lambda (gene 1): 0.675569
#> epsilon (gene 1): -0.362294
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 833.813
#> delta (gene 1): 0.67652
#> phi (cell 1): 0.8723
#> s (cell 1): 0.557243
#> nu (cell 1): 0.404995
#> theta (batch 1): 0.188278
#> betas:    0.8169  -0.0962   0.8079   0.6236   0.7673   0.6540   0.1686   0.2577   0.1113   0.2523  -1.0655   0.2075
#> 
#> sigma: 0.229491
#> lambda (gene 1): 0.753445
#> epsilon (gene 1): -0.201206
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 575.456
#> delta (gene 1): 0.807282
#> phi (cell 1): 1.03556
#> s (cell 1): 0.763932
#> nu (cell 1): 0.407753
#> theta (batch 1): 0.233915
#> betas:    0.6128  -0.0505   0.5951   1.3916   0.1970   0.3139   0.7890   0.2880  -0.3525  -0.1052  -0.2415  -0.4932
#> 
#> sigma: 0.182621
#> lambda (gene 1): 0.762799
#> epsilon (gene 1): 0.081543
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 627.329
#> delta (gene 1): 0.576951
#> phi (cell 1): 0.916332
#> s (cell 1): 0.310997
#> nu (cell 1): 0.427714
#> theta (batch 1): 0.17843
#> betas:    0.9837  -0.0421   0.9221   0.2954   0.9017   0.4305  -0.0790   0.0245  -0.0520  -0.4643  -0.5736  -0.3504
#> 
#> sigma: 0.138402
#> lambda (gene 1): 1.38499
#> epsilon (gene 1): -0.221834
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 702.782
#> delta (gene 1): 0.891559
#> phi (cell 1): 0.989379
#> s (cell 1): 0.473382
#> nu (cell 1): 0.411291
#> theta (batch 1): 0.204428
#> betas:    0.5977   0.0060   1.1726   0.6339   0.2744   0.7348  -0.0415   0.3181   0.0454  -0.5341  -0.6064  -0.2599
#> 
#> sigma: 0.120889
#> lambda (gene 1): 1.60378
#> epsilon (gene 1): 0.220104
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 471.215
#> delta (gene 1): 0.950727
#> phi (cell 1): 0.945551
#> s (cell 1): 0.371723
#> nu (cell 1): 0.408419
#> theta (batch 1): 0.297702
#> betas:    1.0008  -0.1418   0.5184   0.8506   0.8814   0.3026   0.1514   0.2309   0.1172   0.0799  -0.5638   0.0707
#> 
#> sigma: 0.128397
#> lambda (gene 1): 0.549929
#> epsilon (gene 1): 0.0912961
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 671.701
#> delta (gene 1): 0.502086
#> phi (cell 1): 0.918199
#> s (cell 1): 1.56231
#> nu (cell 1): 0.427825
#> theta (batch 1): 0.2961
#> betas:    1.4295  -0.2036   0.5682   0.7487   0.3414   0.2294   0.6450  -0.0379   0.1325  -0.1918  -0.3151   0.2411
#> 
#> sigma: 0.146836
#> lambda (gene 1): 0.392562
#> epsilon (gene 1): -0.503573
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 775.112
#> delta (gene 1): 0.540318
#> phi (cell 1): 0.87002
#> s (cell 1): 0.621966
#> nu (cell 1): 0.426395
#> theta (batch 1): 0.196708
#> betas:    0.3182  -0.0822   0.8373   1.1161   0.5057   0.9858   0.2688   0.3443   0.1134   0.1907  -0.5869   0.4239
#> 
#> sigma: 0.192279
#> lambda (gene 1): 0.918897
#> epsilon (gene 1): -0.411043
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 622.384
#> delta (gene 1): 0.449841
#> phi (cell 1): 0.895026
#> s (cell 1): 0.832375
#> nu (cell 1): 0.419124
#> theta (batch 1): 0.38166
#> betas:    0.3713   0.0925   1.1828   0.7170   0.3571   0.7757   0.0434   0.3074  -0.2567  -0.3820  -0.8561  -0.7971
#> 
#> sigma: 0.15971
#> lambda (gene 1): 0.556781
#> epsilon (gene 1): -0.422263
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 614.318
#> delta (gene 1): 0.684828
#> phi (cell 1): 0.908478
#> s (cell 1): 0.692358
#> nu (cell 1): 0.424082
#> theta (batch 1): 0.0969787
#> betas:    0.9187  -0.0798   0.7934   0.4964   0.5123   0.5146   0.3865   0.0031  -0.2356  -0.1586  -0.5105  -0.0829
#> 
#> sigma: 0.134376
#> lambda (gene 1): 1.21524
#> epsilon (gene 1): -0.0573227
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 700.663
#> delta (gene 1): 0.726724
#> phi (cell 1): 0.773186
#> s (cell 1): 0.47966
#> nu (cell 1): 0.418358
#> theta (batch 1): 0.0571528
#> betas:    0.5656   0.0066   0.9942   0.9356   0.7586  -0.0435   0.4600  -0.0861   0.0255  -0.6655  -0.0743  -0.5553
#> 
#> sigma: 0.213214
#> lambda (gene 1): 0.593614
#> epsilon (gene 1): -0.0615317
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 763.554
#> delta (gene 1): 0.738919
#> phi (cell 1): 0.730294
#> s (cell 1): 0.339294
#> nu (cell 1): 0.429018
#> theta (batch 1): 0.134471
#> betas:    1.1539  -0.1191   0.6499   0.9035   0.1151   0.4780   0.0740   0.3048   0.0873  -0.6240  -0.0613  -0.4984
#> 
#> sigma: 0.137655
#> lambda (gene 1): 2.40557
#> epsilon (gene 1): 0.024339
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 694.039
#> delta (gene 1): 0.486732
#> phi (cell 1): 0.645667
#> s (cell 1): 1.16406
#> nu (cell 1): 0.423002
#> theta (batch 1): 0.329751
#> betas:    0.7729  -0.0757   0.4980   0.6085   0.7539   0.2628   0.7564  -0.1256   0.1595  -0.0565  -0.5876  -0.1707
#> 
#> sigma: 0.138637
#> lambda (gene 1): 0.308092
#> epsilon (gene 1): -0.538666
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 585.109
#> delta (gene 1): 0.579917
#> phi (cell 1): 0.674462
#> s (cell 1): 0.929924
#> nu (cell 1): 0.427282
#> theta (batch 1): 0.333446
#> betas:    0.9884  -0.0569   0.5427   0.5836   0.7602   0.2493   0.6125  -0.0711  -0.3519  -0.0400  -0.7429  -0.1764
#> 
#> sigma: 0.136758
#> lambda (gene 1): 1.36907
#> epsilon (gene 1): -0.275272
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 726.356
#> delta (gene 1): 0.560145
#> phi (cell 1): 0.713524
#> s (cell 1): 0.733439
#> nu (cell 1): 0.417796
#> theta (batch 1): 0.335912
#> betas:    0.0336   0.0506   1.3354   0.7933   0.8507   0.7759   0.3505   0.3128  -0.1891  -0.3619   0.0637  -0.8865
#> 
#> sigma: 0.145369
#> lambda (gene 1): 1.10784
#> epsilon (gene 1): -0.336944
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 966.175
#> delta (gene 1): 0.881255
#> phi (cell 1): 0.698794
#> s (cell 1): 0.310351
#> nu (cell 1): 0.421082
#> theta (batch 1): 0.11671
#> betas:    0.9411  -0.1578   1.0138   0.9132  -0.0656   1.0270   0.6439  -0.3800   0.7907  -0.5417  -0.0883   0.0063
#> 
#> sigma: 0.180495
#> lambda (gene 1): 0.310289
#> epsilon (gene 1): 0.129317
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 759.55
#> delta (gene 1): 0.795983
#> phi (cell 1): 0.75715
#> s (cell 1): 0.398843
#> nu (cell 1): 0.427548
#> theta (batch 1): 0.0307113
#> betas:    0.2747  -0.1274   0.7616   1.1124   0.9434   0.5580   0.4807   0.3648   0.6056  -0.6868   0.1842   0.6163
#> 
#> sigma: 0.194154
#> lambda (gene 1): 1.73581
#> epsilon (gene 1): 0.219985
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 762.223
#> delta (gene 1): 0.794835
#> phi (cell 1): 0.757444
#> s (cell 1): 0.421506
#> nu (cell 1): 0.429476
#> theta (batch 1): 0.0315063
#> betas:    0.9054  -0.0841   0.8864   0.8345   0.4690   0.7115   0.0937  -0.0336   0.2809  -0.2127  -0.5182  -0.1069
#> 
#> sigma: 0.193027
#> lambda (gene 1): 1.84758
#> epsilon (gene 1): -0.117517
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 578.055
#> delta (gene 1): 0.84616
#> phi (cell 1): 0.737303
#> s (cell 1): 0.326206
#> nu (cell 1): 0.41107
#> theta (batch 1): 0.0408439
#> betas:    0.8134  -0.0949   0.4709   0.4828   0.7320   0.7759   0.1605   0.1081   0.0690  -0.0731  -0.2489  -0.3509
#> 
#> sigma: 0.136136
#> lambda (gene 1): 0.995715
#> epsilon (gene 1): -0.0595241
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 909.353
#> delta (gene 1): 0.492885
#> phi (cell 1): 0.711567
#> s (cell 1): 1.00556
#> nu (cell 1): 0.42474
#> theta (batch 1): 0.411137
#> betas:    0.5632  -0.0325   0.9500   1.0749   0.2380   0.6827   0.7306  -0.1360  -0.1787  -0.0495  -0.4357  -0.1599
#> 
#> sigma: 0.158873
#> lambda (gene 1): 0.47682
#> epsilon (gene 1): -0.508774
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 664.735
#> delta (gene 1): 0.612976
#> phi (cell 1): 0.713325
#> s (cell 1): 0.844981
#> nu (cell 1): 0.418323
#> theta (batch 1): 0.145602
#> betas:    0.6717  -0.0885   0.4765   0.6003   1.1110   0.5422   0.4386   0.2596   0.2075  -0.3192  -0.1971  -0.0192
#> 
#> sigma: 0.156552
#> lambda (gene 1): 1.46053
#> epsilon (gene 1): -0.304948
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 634.031
#> delta (gene 1): 0.755813
#> phi (cell 1): 0.741775
#> s (cell 1): 0.36023
#> nu (cell 1): 0.440513
#> theta (batch 1): 0.210132
#> betas:    1.5881  -0.2800   0.4432   0.7282   0.3985   0.1710   0.7502  -0.1131   0.5113  -0.4651  -0.2148   0.4015
#> 
#> sigma: 0.191742
#> lambda (gene 1): 0.708166
#> epsilon (gene 1): 0.170564
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 725.454
#> delta (gene 1): 0.729712
#> phi (cell 1): 0.726157
#> s (cell 1): 0.590296
#> nu (cell 1): 0.419985
#> theta (batch 1): 0.378244
#> betas:    1.0363  -0.0957   0.1529   0.7857   0.5965   0.2036   0.6629   0.0496  -0.1966  -0.2307  -0.3364  -0.3274
#> 
#> sigma: 0.179803
#> lambda (gene 1): 1.0084
#> epsilon (gene 1): -0.0394077
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 738.233
#> delta (gene 1): 0.612687
#> phi (cell 1): 0.759075
#> s (cell 1): 0.908522
#> nu (cell 1): 0.430603
#> theta (batch 1): 0.174589
#> betas:    0.9408  -0.2060   0.3963   0.7494   0.8700   0.6089   0.4119   0.3151   0.2477   0.0707  -0.1456   0.0230
#> 
#> sigma: 0.163647
#> lambda (gene 1): 1.6593
#> epsilon (gene 1): -0.268959
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 663.327
#> delta (gene 1): 0.813941
#> phi (cell 1): 0.728378
#> s (cell 1): 0.476934
#> nu (cell 1): 0.435818
#> theta (batch 1): 0.28807
#> betas:    0.9553  -0.1279   0.8409   0.3196   0.6482   0.5491   0.7197  -0.4244   0.4015  -0.4977  -0.1480  -0.2835
#> 
#> sigma: 0.15474
#> lambda (gene 1): 0.394847
#> epsilon (gene 1): 0.183828
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 864.021
#> delta (gene 1): 1.04346
#> phi (cell 1): 0.690353
#> s (cell 1): 0.577925
#> nu (cell 1): 0.414336
#> theta (batch 1): 0.308718
#> betas:    1.3507  -0.2502   0.5193   0.5240   0.9268   0.1403   0.5092   0.3936  -0.0430  -0.0673   0.0301   0.1543
#> 
#> sigma: 0.170993
#> lambda (gene 1): 1.52484
#> epsilon (gene 1): 0.282732
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 706.943
#> delta (gene 1): 0.691779
#> phi (cell 1): 0.763777
#> s (cell 1): 0.275736
#> nu (cell 1): 0.426857
#> theta (batch 1): 0.396597
#> betas:    0.4466   0.0414   1.0459   1.1351   0.6585   0.5397   0.3412   0.2622  -0.4827  -0.1014  -0.3863  -0.7999
#> 
#> sigma: 0.137286
#> lambda (gene 1): 1.03506
#> epsilon (gene 1): -0.217497
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 688.185
#> delta (gene 1): 0.561592
#> phi (cell 1): 0.839299
#> s (cell 1): 0.564464
#> nu (cell 1): 0.428266
#> theta (batch 1): 0.0871192
#> betas:   -0.0686   0.0062   0.6886   0.9388   0.5979   1.4585   0.0431   0.2025  -0.0063   0.6210  -0.8485  -0.2051
#> 
#> sigma: 0.239362
#> lambda (gene 1): 2.32432
#> epsilon (gene 1): -0.505004
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 718.793
#> delta (gene 1): 0.837255
#> phi (cell 1): 0.920246
#> s (cell 1): 0.219163
#> nu (cell 1): 0.413073
#> theta (batch 1): 0.365542
#> betas:    1.1360   0.0402   0.5255   1.0198   0.1834   0.4516  -0.1307   0.1155  -0.6282  -0.5467  -0.8207  -0.5412
#> 
#> sigma: 0.197814
#> lambda (gene 1): 3.32062
#> epsilon (gene 1): 0.133002
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 16980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 787.272
#> delta (gene 1): 0.573575
#> phi (cell 1): 0.884912
#> s (cell 1): 0.567325
#> nu (cell 1): 0.419762
#> theta (batch 1): 0.0362566
#> betas:    1.1960  -0.1087   1.5947   0.5751   0.2679   0.5942   0.0897   0.2714  -0.0062  -0.5529  -0.0530  -0.5639
#> 
#> sigma: 0.234799
#> lambda (gene 1): 1.99407
#> epsilon (gene 1): -0.3315
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 643.594
#> delta (gene 1): 0.705621
#> phi (cell 1): 0.828158
#> s (cell 1): 0.447798
#> nu (cell 1): 0.418774
#> theta (batch 1): 0.014263
#> betas:    1.0414  -0.1617   0.7415   0.2973   0.8943   0.3014   0.8036   0.0794  -0.0888  -0.0167  -0.3199   0.0476
#> 
#> sigma: 0.173716
#> lambda (gene 1): 1.2324
#> epsilon (gene 1): -0.055031
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 673.494
#> delta (gene 1): 0.821319
#> phi (cell 1): 0.757315
#> s (cell 1): 0.407061
#> nu (cell 1): 0.420454
#> theta (batch 1): 0.0142368
#> betas:    1.4464  -0.1878   0.4670   1.0113   0.1765   0.5228   0.3255   0.0431   0.1089  -0.1067  -0.2721  -0.1983
#> 
#> sigma: 0.156312
#> lambda (gene 1): 0.564714
#> epsilon (gene 1): -0.114309
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 701.545
#> delta (gene 1): 0.659193
#> phi (cell 1): 0.760592
#> s (cell 1): 0.455508
#> nu (cell 1): 0.417088
#> theta (batch 1): 0.010386
#> betas:    0.8565  -0.1272   0.5989   0.3723   1.2803   0.2230   0.5367   0.0206   0.4538  -0.5205  -0.0085  -0.2832
#> 
#> sigma: 0.138597
#> lambda (gene 1): 0.794545
#> epsilon (gene 1): -0.137323
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 662.643
#> delta (gene 1): 0.528835
#> phi (cell 1): 0.813399
#> s (cell 1): 0.4947
#> nu (cell 1): 0.410519
#> theta (batch 1): 0.055765
#> betas:   -0.2249   0.0347   0.0251   1.7130   0.6659   0.9062   0.2912   0.2714   0.3025  -0.1531  -0.4197  -0.2777
#> 
#> sigma: 0.155122
#> lambda (gene 1): 1.07292
#> epsilon (gene 1): -0.355508
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 616.678
#> delta (gene 1): 0.728631
#> phi (cell 1): 0.76655
#> s (cell 1): 0.581884
#> nu (cell 1): 0.421286
#> theta (batch 1): 0.0822276
#> betas:    0.8661  -0.1003   0.3948   0.6671   0.5146   0.9488   0.2663  -0.1113   0.1086  -0.1110  -0.3390  -0.2675
#> 
#> sigma: 0.114277
#> lambda (gene 1): 1.52128
#> epsilon (gene 1): -0.127089
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 661.198
#> delta (gene 1): 0.447528
#> phi (cell 1): 0.806687
#> s (cell 1): 0.72251
#> nu (cell 1): 0.407646
#> theta (batch 1): 0.0690069
#> betas:    0.2422   0.0927   0.6408   1.0451   0.7708   0.7326   0.0270   0.2636  -0.3416   0.0055  -0.9560  -0.7125
#> 
#> sigma: 0.215142
#> lambda (gene 1): 0.8053
#> epsilon (gene 1): -0.568619
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 891.734
#> delta (gene 1): 0.968984
#> phi (cell 1): 0.722463
#> s (cell 1): 0.487728
#> nu (cell 1): 0.416785
#> theta (batch 1): 0.0863727
#> betas:    0.7293  -0.0281   0.8102   0.7314   0.6217   0.5690   0.1390   0.3786  -0.0994  -0.5857  -0.3253  -0.2960
#> 
#> sigma: 0.126901
#> lambda (gene 1): 0.0843225
#> epsilon (gene 1): 0.322029
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 633.63
#> delta (gene 1): 0.856023
#> phi (cell 1): 0.759325
#> s (cell 1): 0.705898
#> nu (cell 1): 0.426246
#> theta (batch 1): 0.406173
#> betas:    0.5049  -0.0297   1.0673   0.4546   0.5399   0.9799   0.0223   0.0612   0.3383  -0.1563  -0.3638  -0.3728
#> 
#> sigma: 0.16707
#> lambda (gene 1): 0.642003
#> epsilon (gene 1): -0.156646
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 806.293
#> delta (gene 1): 0.65279
#> phi (cell 1): 0.898227
#> s (cell 1): 0.488709
#> nu (cell 1): 0.415774
#> theta (batch 1): 0.397717
#> betas:    1.1300  -0.1767   0.9694   0.6763   0.4835   0.9186   0.0351   0.1346   0.3221  -0.0750  -0.5087   0.2339
#> 
#> sigma: 0.171062
#> lambda (gene 1): 0.749227
#> epsilon (gene 1): -0.227457
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 798.983
#> delta (gene 1): 0.643559
#> phi (cell 1): 0.876388
#> s (cell 1): 0.383779
#> nu (cell 1): 0.416439
#> theta (batch 1): 0.153505
#> betas:    0.6505   0.0215   0.7365   1.0139   0.2897   0.6223  -0.0276   0.1566   0.0725  -0.5556  -0.4420  -0.4439
#> 
#> sigma: 0.123508
#> lambda (gene 1): 0.597848
#> epsilon (gene 1): -0.29602
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 874.302
#> delta (gene 1): 0.840647
#> phi (cell 1): 0.973244
#> s (cell 1): 0.731902
#> nu (cell 1): 0.419177
#> theta (batch 1): 0.217144
#> betas:    0.0950   0.0070   0.6660   0.8932   0.9265   0.5292   0.5698   0.3753  -0.2210   0.0329  -0.1366  -0.6411
#> 
#> sigma: 0.140665
#> lambda (gene 1): 1.10163
#> epsilon (gene 1): -0.0134478
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 678.072
#> delta (gene 1): 0.525047
#> phi (cell 1): 1.01238
#> s (cell 1): 1.23288
#> nu (cell 1): 0.428184
#> theta (batch 1): 0.334116
#> betas:    0.7396  -0.0940   0.3894   0.5098   0.9622   0.5650   0.2890   0.4862  -0.0963  -0.0885  -0.5521   0.3643
#> 
#> sigma: 0.153115
#> lambda (gene 1): 1.09439
#> epsilon (gene 1): -0.407647
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 621.75
#> delta (gene 1): 1.01906
#> phi (cell 1): 0.907554
#> s (cell 1): 0.628107
#> nu (cell 1): 0.412447
#> theta (batch 1): 0.212805
#> betas:    1.2759  -0.1385   0.6744   0.5681   0.6726   0.4136   0.2168   0.0104   0.0518  -0.0198  -0.7643   0.2780
#> 
#> sigma: 0.196941
#> lambda (gene 1): 0.523403
#> epsilon (gene 1): 0.141517
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 827.435
#> delta (gene 1): 0.789551
#> phi (cell 1): 0.757579
#> s (cell 1): 0.365405
#> nu (cell 1): 0.443936
#> theta (batch 1): 0.362069
#> betas:    0.6825  -0.0657   0.8749   0.7243   0.6227   0.7216   0.3265   0.0753   0.4585  -0.5630  -0.4267   0.0025
#> 
#> sigma: 0.12738
#> lambda (gene 1): 0.530764
#> epsilon (gene 1): 0.0771575
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 589.605
#> delta (gene 1): 0.730935
#> phi (cell 1): 0.711627
#> s (cell 1): 1.30304
#> nu (cell 1): 0.410396
#> theta (batch 1): 0.421201
#> betas:    0.8482  -0.0539   0.4892   0.4743   0.9806   0.9026  -0.5482   0.5918   0.2370  -0.6934  -0.5510  -0.0107
#> 
#> sigma: 0.169752
#> lambda (gene 1): 0.855426
#> epsilon (gene 1): 0.0483365
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 783.608
#> delta (gene 1): 0.526683
#> phi (cell 1): 0.607873
#> s (cell 1): 0.472817
#> nu (cell 1): 0.430773
#> theta (batch 1): 0.23895
#> betas:    0.9456  -0.0828   0.8588   0.6038   0.4791   0.6862   0.3052  -0.2164   0.0280  -0.0758  -0.4856  -0.4492
#> 
#> sigma: 0.134717
#> lambda (gene 1): 1.01049
#> epsilon (gene 1): -0.443668
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 771.942
#> delta (gene 1): 0.564211
#> phi (cell 1): 0.643452
#> s (cell 1): 0.558008
#> nu (cell 1): 0.423666
#> theta (batch 1): 0.283744
#> betas:    1.1780  -0.0488   0.1904   1.0854   0.2897   0.8517  -0.1526  -0.2416   0.1206  -0.4235  -0.8174  -0.4033
#> 
#> sigma: 0.268006
#> lambda (gene 1): 0.409269
#> epsilon (gene 1): -0.287462
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 691.967
#> delta (gene 1): 0.801164
#> phi (cell 1): 0.657872
#> s (cell 1): 0.69494
#> nu (cell 1): 0.428743
#> theta (batch 1): 0.0874131
#> betas:    0.8963  -0.1067   0.3873   0.5595   0.5731   0.9550   0.1649   0.1314   0.0819  -0.3596  -0.2261  -0.2542
#> 
#> sigma: 0.167226
#> lambda (gene 1): 1.86959
#> epsilon (gene 1): 0.107106
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 733.456
#> delta (gene 1): 0.759012
#> phi (cell 1): 0.68366
#> s (cell 1): 0.587752
#> nu (cell 1): 0.410696
#> theta (batch 1): 0.181573
#> betas:    0.8404   0.0271   0.4605   1.0153   0.2375   0.5858   0.1819  -0.0171  -0.5237  -0.0558  -0.7905  -0.9733
#> 
#> sigma: 0.146766
#> lambda (gene 1): 0.819977
#> epsilon (gene 1): -0.0189704
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 923.678
#> delta (gene 1): 0.604111
#> phi (cell 1): 0.673244
#> s (cell 1): 0.273217
#> nu (cell 1): 0.403423
#> theta (batch 1): 0.247591
#> betas:    0.5565  -0.1576   0.6851   1.2636   0.0923   0.9433   0.9755   0.1235  -0.1413   0.3694   0.1689  -0.1637
#> 
#> sigma: 0.184898
#> lambda (gene 1): 0.315448
#> epsilon (gene 1): -0.431257
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 643.7
#> delta (gene 1): 0.523337
#> phi (cell 1): 0.754965
#> s (cell 1): 0.256221
#> nu (cell 1): 0.414566
#> theta (batch 1): 0.170002
#> betas:    0.0762   0.0973   1.4261   0.4411   0.8152   1.0125  -0.3885   0.5595  -0.2794  -0.1134  -0.7077  -0.6090
#> 
#> sigma: 0.165041
#> lambda (gene 1): 0.72457
#> epsilon (gene 1): -0.467305
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 594.359
#> delta (gene 1): 0.555043
#> phi (cell 1): 0.81018
#> s (cell 1): 2.32185
#> nu (cell 1): 0.420642
#> theta (batch 1): 0.252005
#> betas:    1.2197  -0.1844  -0.0075   1.4958   0.1190   0.3262   0.5757   0.2987  -0.2915   0.2872  -0.6970   0.5508
#> 
#> sigma: 0.201017
#> lambda (gene 1): 1.09715
#> epsilon (gene 1): -0.390853
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 794.249
#> delta (gene 1): 0.449752
#> phi (cell 1): 0.730792
#> s (cell 1): 1.32878
#> nu (cell 1): 0.432598
#> theta (batch 1): 0.265331
#> betas:    0.6162  -0.0928  -0.0989   1.3572   0.3240   1.1513   0.0654  -0.0393   0.4729  -0.4058   0.0235  -0.2883
#> 
#> sigma: 0.155829
#> lambda (gene 1): 1.16086
#> epsilon (gene 1): -0.631892
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 566.541
#> delta (gene 1): 0.70937
#> phi (cell 1): 0.760939
#> s (cell 1): 0.510741
#> nu (cell 1): 0.415549
#> theta (batch 1): 0.0950893
#> betas:    0.7464  -0.1743   0.7852   0.4676   0.9209   1.0558   0.0815   0.4526   0.1170   0.2398  -0.1772  -0.0811
#> 
#> sigma: 0.188339
#> lambda (gene 1): 0.86067
#> epsilon (gene 1): -0.187092
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 719.604
#> delta (gene 1): 0.939273
#> phi (cell 1): 0.777064
#> s (cell 1): 0.451265
#> nu (cell 1): 0.418502
#> theta (batch 1): 0.0925969
#> betas:    1.2452  -0.1268   0.3898   0.7336   0.7631   0.1927   0.1750   0.4456  -0.4419  -0.0424  -0.5599  -0.1167
#> 
#> sigma: 0.151072
#> lambda (gene 1): 2.84596
#> epsilon (gene 1): 0.225422
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 617.861
#> delta (gene 1): 0.560637
#> phi (cell 1): 0.856258
#> s (cell 1): 1.0136
#> nu (cell 1): 0.429469
#> theta (batch 1): 0.281915
#> betas:    1.3439  -0.2254  -0.2148   1.0933   0.3333   0.9057  -0.2299   0.5656   0.1869  -0.0220  -0.3423   0.1040
#> 
#> sigma: 0.190864
#> lambda (gene 1): 0.909485
#> epsilon (gene 1): -0.40579
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 949.94
#> delta (gene 1): 0.708387
#> phi (cell 1): 0.730864
#> s (cell 1): 0.444143
#> nu (cell 1): 0.433693
#> theta (batch 1): 0.0101699
#> betas:    0.0524   0.1061   0.5615   1.3340   0.6237   0.6630   0.0655   0.4552  -0.3570  -0.0874  -0.5898  -0.8444
#> 
#> sigma: 0.133326
#> lambda (gene 1): 0.634894
#> epsilon (gene 1): -0.189002
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 580.322
#> delta (gene 1): 0.897272
#> phi (cell 1): 0.674899
#> s (cell 1): 0.418553
#> nu (cell 1): 0.440232
#> theta (batch 1): 0.020737
#> betas:    0.4061  -0.0238   0.6910   1.1645   0.5345   0.5209   0.3064   0.6868  -0.5076   0.1476  -0.5451  -0.2015
#> 
#> sigma: 0.191942
#> lambda (gene 1): 0.568856
#> epsilon (gene 1): 0.143665
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 685.259
#> delta (gene 1): 0.555022
#> phi (cell 1): 0.691652
#> s (cell 1): 0.556748
#> nu (cell 1): 0.429924
#> theta (batch 1): 0.117061
#> betas:    0.6291  -0.0737   0.9893   1.0073   0.6401   0.1582   0.8690   0.2310  -0.3164   0.2738  -0.8245   0.1361
#> 
#> sigma: 0.181997
#> lambda (gene 1): 1.83528
#> epsilon (gene 1): -0.284838
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 713.712
#> delta (gene 1): 0.63419
#> phi (cell 1): 0.811655
#> s (cell 1): 0.329592
#> nu (cell 1): 0.408024
#> theta (batch 1): 0.156563
#> betas:    1.0006  -0.1981   1.0611   0.3062   1.3275   0.0666   0.7037   0.0147   0.2434  -0.0792  -0.0107  -0.0723
#> 
#> sigma: 0.175664
#> lambda (gene 1): 1.50115
#> epsilon (gene 1): -0.215416
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 695.142
#> delta (gene 1): 0.607365
#> phi (cell 1): 0.754863
#> s (cell 1): 0.535058
#> nu (cell 1): 0.412754
#> theta (batch 1): 0.270884
#> betas:    0.5192   0.0089   0.3419   0.8368   0.6922   0.9225  -0.2606   0.3758  -0.1422  -0.2832  -0.4500  -0.3318
#> 
#> sigma: 0.145407
#> lambda (gene 1): 0.580127
#> epsilon (gene 1): -0.335832
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 755.658
#> delta (gene 1): 1.00532
#> phi (cell 1): 0.767671
#> s (cell 1): 0.190094
#> nu (cell 1): 0.423165
#> theta (batch 1): 0.251471
#> betas:    0.7038  -0.0744   0.4617   1.1135   0.6181   0.3704   0.5783   0.1240   0.0564  -0.0591  -0.7294   0.3418
#> 
#> sigma: 0.156575
#> lambda (gene 1): 1.42082
#> epsilon (gene 1): 0.228231
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 691.31
#> delta (gene 1): 0.474723
#> phi (cell 1): 0.726911
#> s (cell 1): 0.516872
#> nu (cell 1): 0.414683
#> theta (batch 1): 0.0953563
#> betas:    1.1757  -0.1909   0.3919   1.2095   0.3736   0.3362   0.5264   0.3696   0.0730  -0.2908   0.0661  -0.1831
#> 
#> sigma: 0.149825
#> lambda (gene 1): 2.30873
#> epsilon (gene 1): -0.502333
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 820.746
#> delta (gene 1): 0.521536
#> phi (cell 1): 0.718551
#> s (cell 1): 1.19245
#> nu (cell 1): 0.42566
#> theta (batch 1): 0.253539
#> betas:    0.5924  -0.0863   0.9594   1.0277   0.7365   0.3395   0.3403   0.2274   0.4261  -0.7692   0.4281  -0.6752
#> 
#> sigma: 0.182076
#> lambda (gene 1): 0.569053
#> epsilon (gene 1): -0.370362
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 704.469
#> delta (gene 1): 0.733457
#> phi (cell 1): 0.705878
#> s (cell 1): 0.563754
#> nu (cell 1): 0.416117
#> theta (batch 1): 0.212441
#> betas:    1.2625  -0.2925   0.4065   1.0487   0.3415   0.6925   0.6141   0.2187   0.2841  -0.0983   0.2648   0.4072
#> 
#> sigma: 0.16338
#> lambda (gene 1): 0.493159
#> epsilon (gene 1): -0.103033
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 734.513
#> delta (gene 1): 0.698212
#> phi (cell 1): 0.770951
#> s (cell 1): 0.307379
#> nu (cell 1): 0.417696
#> theta (batch 1): 0.143743
#> betas:    0.6832  -0.1023   0.1673   1.1347   0.8769   0.6830  -0.2312   0.6131  -0.0297  -0.1808  -0.1937  -0.1843
#> 
#> sigma: 0.175485
#> lambda (gene 1): 1.58419
#> epsilon (gene 1): -0.0847052
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 626.611
#> delta (gene 1): 0.612823
#> phi (cell 1): 0.942235
#> s (cell 1): 0.635662
#> nu (cell 1): 0.424936
#> theta (batch 1): 0.135899
#> betas:    1.1190  -0.0547   0.5212   0.6340   0.7735   0.5463  -0.0642  -0.0544   0.1700  -0.2329  -0.8837  -0.2969
#> 
#> sigma: 0.13582
#> lambda (gene 1): 1.36728
#> epsilon (gene 1): -0.350391
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 607.823
#> delta (gene 1): 0.731895
#> phi (cell 1): 0.924664
#> s (cell 1): 0.320254
#> nu (cell 1): 0.439281
#> theta (batch 1): 0.103976
#> betas:    0.2197   0.0134   1.1842   0.6810   0.9067   0.5503   0.3904   0.2267   0.2231  -0.2571  -0.5600  -0.3614
#> 
#> sigma: 0.156168
#> lambda (gene 1): 0.691123
#> epsilon (gene 1): -0.0400565
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 794.095
#> delta (gene 1): 0.881047
#> phi (cell 1): 0.911547
#> s (cell 1): 0.410142
#> nu (cell 1): 0.419073
#> theta (batch 1): 0.224048
#> betas:    0.5331  -0.2300   1.1877   0.7757   0.8557   0.6423   0.7467   0.1344   0.6934   0.1672   0.0367   0.6652
#> 
#> sigma: 0.186015
#> lambda (gene 1): 2.0834
#> epsilon (gene 1): -0.0235135
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 790.847
#> delta (gene 1): 0.962419
#> phi (cell 1): 0.848554
#> s (cell 1): 0.762009
#> nu (cell 1): 0.402733
#> theta (batch 1): 0.223735
#> betas:    0.9784  -0.0202   1.0745   0.4946   0.6216   0.3016   0.2837   0.1747  -0.1959  -0.1593  -0.8750  -0.2647
#> 
#> sigma: 0.186821
#> lambda (gene 1): 0.610389
#> epsilon (gene 1): 0.100567
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 690.533
#> delta (gene 1): 1.10795
#> phi (cell 1): 0.788126
#> s (cell 1): 0.345832
#> nu (cell 1): 0.407493
#> theta (batch 1): 0.196031
#> betas:    0.3964   0.0097   0.8379   1.3647   0.1172   1.0538  -0.1243   0.1889  -0.1415   0.1587  -0.6803  -0.2458
#> 
#> sigma: 0.18504
#> lambda (gene 1): 1.50013
#> epsilon (gene 1): 0.127741
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 652.694
#> delta (gene 1): 0.596
#> phi (cell 1): 0.685867
#> s (cell 1): 0.422494
#> nu (cell 1): 0.409279
#> theta (batch 1): 0.11794
#> betas:    1.1569  -0.1509   0.5986   0.3837   0.3620   0.9365  -0.0422   0.3362   0.1397  -0.2102  -0.3363  -0.2037
#> 
#> sigma: 0.161101
#> lambda (gene 1): 0.566959
#> epsilon (gene 1): -0.307833
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 741.011
#> delta (gene 1): 0.668572
#> phi (cell 1): 0.758792
#> s (cell 1): 0.153314
#> nu (cell 1): 0.422309
#> theta (batch 1): 0.257696
#> betas:    0.8768  -0.0681   0.4163   0.9602   0.4235   0.4771   0.3627   0.0182   0.2921  -0.7257  -0.2193  -0.0357
#> 
#> sigma: 0.128529
#> lambda (gene 1): 0.988408
#> epsilon (gene 1): -0.133758
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 540.432
#> delta (gene 1): 0.637261
#> phi (cell 1): 0.74602
#> s (cell 1): 0.633118
#> nu (cell 1): 0.430269
#> theta (batch 1): 0.373947
#> betas:    1.2447  -0.0637   1.0712   0.3382   0.8530  -0.0212   0.2249  -0.1496   0.0585  -0.4960  -0.6065  -0.2532
#> 
#> sigma: 0.186202
#> lambda (gene 1): 0.451533
#> epsilon (gene 1): -0.255186
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 708.159
#> delta (gene 1): 0.48503
#> phi (cell 1): 0.660886
#> s (cell 1): 0.512582
#> nu (cell 1): 0.438215
#> theta (batch 1): 0.228475
#> betas:    0.4909  -0.0658   0.2020   0.8586   0.4067   0.8952   0.4475   0.1972   0.0981  -0.2207  -0.3859   0.1382
#> 
#> sigma: 0.118568
#> lambda (gene 1): 0.643083
#> epsilon (gene 1): -0.42644
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 797.448
#> delta (gene 1): 0.766872
#> phi (cell 1): 0.743489
#> s (cell 1): 0.309273
#> nu (cell 1): 0.421418
#> theta (batch 1): 0.301444
#> betas:    0.5672  -0.0351   1.3156   0.8467   0.8573   0.1002   0.4201   0.1843   0.1822  -0.4402  -0.3470  -0.0828
#> 
#> sigma: 0.137536
#> lambda (gene 1): 0.996251
#> epsilon (gene 1): -0.0440467
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 705.019
#> delta (gene 1): 0.606867
#> phi (cell 1): 0.898602
#> s (cell 1): 0.271997
#> nu (cell 1): 0.420804
#> theta (batch 1): 0.272818
#> betas:    0.3193   0.0102   0.1282   0.8607   1.2941   0.2786   0.3718  -0.0689   0.1743  -0.0521  -0.6506  -0.2379
#> 
#> sigma: 0.151214
#> lambda (gene 1): 0.772662
#> epsilon (gene 1): -0.381456
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 712.351
#> delta (gene 1): 1.13404
#> phi (cell 1): 0.990439
#> s (cell 1): 0.321202
#> nu (cell 1): 0.418804
#> theta (batch 1): 0.256146
#> betas:    0.4263   0.0011   0.8286   0.7369   0.8057   0.6464   0.1830  -0.0951   0.0761  -0.1158  -0.8038  -0.2539
#> 
#> sigma: 0.149791
#> lambda (gene 1): 1.36762
#> epsilon (gene 1): 0.457421
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 640.537
#> delta (gene 1): 0.860658
#> phi (cell 1): 1.00996
#> s (cell 1): 0.552329
#> nu (cell 1): 0.412906
#> theta (batch 1): 0.0943342
#> betas:    0.3316   0.0795   0.1128   1.0286   0.5577   0.7393  -0.0406   0.1059  -0.0837  -0.3426  -0.8784  -0.3493
#> 
#> sigma: 0.117021
#> lambda (gene 1): 0.66056
#> epsilon (gene 1): 0.133596
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 17980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 622.032
#> delta (gene 1): 0.819528
#> phi (cell 1): 1.04794
#> s (cell 1): 0.283247
#> nu (cell 1): 0.420004
#> theta (batch 1): 0.259229
#> betas:    0.5847  -0.0047   0.6047   1.0592   0.4687   0.7821   0.1342   0.0692  -0.0481  -0.2623  -0.8075  -0.0326
#> 
#> sigma: 0.182967
#> lambda (gene 1): 1.36739
#> epsilon (gene 1): 0.143449
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 705.571
#> delta (gene 1): 0.573652
#> phi (cell 1): 0.951649
#> s (cell 1): 0.565921
#> nu (cell 1): 0.411024
#> theta (batch 1): 0.156763
#> betas:    0.7164  -0.1054   0.5213   1.0316   0.7777   0.2109   0.3822   0.5583   0.0113  -0.4505   0.1948  -0.1905
#> 
#> sigma: 0.15802
#> lambda (gene 1): 0.460005
#> epsilon (gene 1): -0.354585
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 663.255
#> delta (gene 1): 0.756543
#> phi (cell 1): 1.0078
#> s (cell 1): 0.158774
#> nu (cell 1): 0.414114
#> theta (batch 1): 0.333918
#> betas:    0.9029  -0.1335   0.2656   0.9355   0.2908   0.5495   1.0422  -0.2642   0.2873  -0.3455   0.0112  -0.2198
#> 
#> sigma: 0.176055
#> lambda (gene 1): 1.91395
#> epsilon (gene 1): -0.0619639
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 562.639
#> delta (gene 1): 0.65185
#> phi (cell 1): 0.863333
#> s (cell 1): 0.567802
#> nu (cell 1): 0.415937
#> theta (batch 1): 0.677691
#> betas:    0.6277  -0.0547   0.6544   1.1678   0.6589   0.3610   0.5634  -0.0441  -0.0038   0.1155  -0.3112  -0.6893
#> 
#> sigma: 0.163262
#> lambda (gene 1): 0.623956
#> epsilon (gene 1): -0.376962
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 666.252
#> delta (gene 1): 0.53684
#> phi (cell 1): 0.889165
#> s (cell 1): 0.394237
#> nu (cell 1): 0.415985
#> theta (batch 1): 0.429984
#> betas:    0.4698   0.0346   0.5745   1.2101   0.0421   0.8766   0.0022   0.4727  -0.4712  -0.1471  -0.7094  -0.4894
#> 
#> sigma: 0.158576
#> lambda (gene 1): 0.741809
#> epsilon (gene 1): -0.288701
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 585.217
#> delta (gene 1): 0.609334
#> phi (cell 1): 0.822177
#> s (cell 1): 0.732137
#> nu (cell 1): 0.417305
#> theta (batch 1): 0.0591331
#> betas:    0.6262  -0.1017   1.0468   0.8957   0.5990   0.7112   0.2302   0.5322   0.0137  -0.3353  -0.0269   0.0310
#> 
#> sigma: 0.161806
#> lambda (gene 1): 0.63221
#> epsilon (gene 1): -0.230684
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 615.702
#> delta (gene 1): 0.601768
#> phi (cell 1): 0.751887
#> s (cell 1): 0.462811
#> nu (cell 1): 0.425115
#> theta (batch 1): 0.128004
#> betas:    1.2691  -0.1099   1.0516   0.7881  -0.2838   0.7339   0.4609  -0.3597   0.2178  -0.3753  -0.8546   0.0090
#> 
#> sigma: 0.137748
#> lambda (gene 1): 0.738635
#> epsilon (gene 1): -0.106423
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 648.809
#> delta (gene 1): 0.814426
#> phi (cell 1): 0.768157
#> s (cell 1): 0.771479
#> nu (cell 1): 0.419712
#> theta (batch 1): 0.242849
#> betas:    0.7851  -0.0558  -0.0625   1.0062   0.8491   0.2329   0.1198   0.4140  -0.1628  -0.5640   0.2278  -0.5867
#> 
#> sigma: 0.150911
#> lambda (gene 1): 1.76023
#> epsilon (gene 1): -0.0471134
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 767.213
#> delta (gene 1): 0.789473
#> phi (cell 1): 0.907756
#> s (cell 1): 0.306317
#> nu (cell 1): 0.423288
#> theta (batch 1): 0.532068
#> betas:    0.7254  -0.0117   0.5282   0.7406   0.7315   0.6215   0.1171   0.1069  -0.0343  -0.4797  -0.5733  -0.1797
#> 
#> sigma: 0.165528
#> lambda (gene 1): 0.589735
#> epsilon (gene 1): 0.0767064
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 803.637
#> delta (gene 1): 0.794737
#> phi (cell 1): 0.783685
#> s (cell 1): 0.35857
#> nu (cell 1): 0.425175
#> theta (batch 1): 0.308104
#> betas:    0.7722  -0.1500   0.2421   1.4651   0.2604   0.7838   0.3844   0.2384   0.1291  -0.1288   0.1869  -0.1890
#> 
#> sigma: 0.186534
#> lambda (gene 1): 0.939051
#> epsilon (gene 1): -0.112411
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 701.914
#> delta (gene 1): 1.43285
#> phi (cell 1): 0.84257
#> s (cell 1): 0.594209
#> nu (cell 1): 0.427453
#> theta (batch 1): 0.319005
#> betas:    0.8061   0.0252   0.5435   0.6358   0.7547   0.5568   0.0766  -0.0087   0.0198  -0.6405  -0.4277  -0.5640
#> 
#> sigma: 0.162964
#> lambda (gene 1): 0.676209
#> epsilon (gene 1): 0.49852
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 754.101
#> delta (gene 1): 0.83976
#> phi (cell 1): 0.906753
#> s (cell 1): 0.570567
#> nu (cell 1): 0.408499
#> theta (batch 1): 0.118781
#> betas:    1.0540  -0.0430   1.0472   0.3239   0.3530   0.6052   0.3449  -0.3193  -0.0688  -0.2446  -0.3941  -0.7476
#> 
#> sigma: 0.139295
#> lambda (gene 1): 1.72566
#> epsilon (gene 1): -0.0819375
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 862.952
#> delta (gene 1): 0.611745
#> phi (cell 1): 0.921412
#> s (cell 1): 0.438306
#> nu (cell 1): 0.412733
#> theta (batch 1): 0.060867
#> betas:    1.0054  -0.1805   0.5316   0.5841   1.0654   0.5426   0.5276   0.2022   0.0503   0.0037  -0.0739  -0.0699
#> 
#> sigma: 0.164033
#> lambda (gene 1): 0.555205
#> epsilon (gene 1): -0.317723
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 832.366
#> delta (gene 1): 0.814407
#> phi (cell 1): 0.907567
#> s (cell 1): 0.53952
#> nu (cell 1): 0.418401
#> theta (batch 1): 0.123456
#> betas:    0.3805   0.0222   0.9046   0.4598   0.9030   0.6897   0.0118   0.3725  -0.1157  -0.0931  -0.9281   0.2251
#> 
#> sigma: 0.186681
#> lambda (gene 1): 0.827063
#> epsilon (gene 1): -0.00468452
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 699.477
#> delta (gene 1): 0.709907
#> phi (cell 1): 0.765304
#> s (cell 1): 0.484943
#> nu (cell 1): 0.403185
#> theta (batch 1): 0.121351
#> betas:    0.2025   0.0545   1.1907   0.6981   0.7152   0.4008   0.5512  -0.1788   0.3538  -0.4908  -0.5307  -0.2814
#> 
#> sigma: 0.239293
#> lambda (gene 1): 1.46071
#> epsilon (gene 1): -0.144501
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 766.982
#> delta (gene 1): 0.672915
#> phi (cell 1): 0.716861
#> s (cell 1): 0.393911
#> nu (cell 1): 0.420164
#> theta (batch 1): 0.045838
#> betas:    1.0535  -0.1070   0.9153   0.5990   0.6130   0.3647   0.4549   0.1269  -0.0060  -0.1768  -0.0773  -0.9212
#> 
#> sigma: 0.13503
#> lambda (gene 1): 2.53746
#> epsilon (gene 1): -0.273235
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 520.624
#> delta (gene 1): 0.700619
#> phi (cell 1): 0.726308
#> s (cell 1): 0.351626
#> nu (cell 1): 0.427193
#> theta (batch 1): 0.0896628
#> betas:    1.3363  -0.1505   0.3550   0.6136   0.4434   0.6832  -0.3204   0.7210  -0.2150  -0.3034  -0.1687  -0.4311
#> 
#> sigma: 0.113844
#> lambda (gene 1): 1.3392
#> epsilon (gene 1): -0.137803
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 743.045
#> delta (gene 1): 0.493686
#> phi (cell 1): 0.648931
#> s (cell 1): 0.606096
#> nu (cell 1): 0.422156
#> theta (batch 1): 0.651794
#> betas:    0.8571  -0.0711   0.7385   0.8759   0.0519   1.1751  -0.1038  -0.1244   0.2786  -0.0884  -0.7769  -0.1217
#> 
#> sigma: 0.134158
#> lambda (gene 1): 1.02373
#> epsilon (gene 1): -0.528021
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 650.408
#> delta (gene 1): 0.542473
#> phi (cell 1): 0.663462
#> s (cell 1): 0.361454
#> nu (cell 1): 0.415959
#> theta (batch 1): 0.790149
#> betas:    0.2618  -0.0911   0.3722   1.3679   1.0175   0.1025   0.9899   0.1842   0.5032  -0.4572   0.4801  -0.3101
#> 
#> sigma: 0.172738
#> lambda (gene 1): 0.616107
#> epsilon (gene 1): -0.487328
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 805.59
#> delta (gene 1): 0.794491
#> phi (cell 1): 0.595425
#> s (cell 1): 0.204464
#> nu (cell 1): 0.420512
#> theta (batch 1): 0.3669
#> betas:    1.0497  -0.0936   0.3795   0.6975   0.4028   0.4995   0.5388  -0.2451   0.1012  -0.3381  -0.6784   0.1597
#> 
#> sigma: 0.165758
#> lambda (gene 1): 2.48712
#> epsilon (gene 1): 0.116376
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 690.294
#> delta (gene 1): 0.944759
#> phi (cell 1): 0.6589
#> s (cell 1): 0.676792
#> nu (cell 1): 0.415834
#> theta (batch 1): 0.349994
#> betas:    0.5875  -0.0896   0.4300   0.7357   0.8934   0.5253   0.2123   0.1198   0.3861  -0.2764  -0.1923  -0.0310
#> 
#> sigma: 0.0994366
#> lambda (gene 1): 0.937957
#> epsilon (gene 1): 0.0983468
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 800.436
#> delta (gene 1): 0.460387
#> phi (cell 1): 0.721012
#> s (cell 1): 0.481539
#> nu (cell 1): 0.432818
#> theta (batch 1): 0.192173
#> betas:    0.7628   0.0154   0.7201   0.6516   0.2881   0.6436   0.3146  -0.3561  -0.1683  -0.3016  -0.6888  -0.3964
#> 
#> sigma: 0.22342
#> lambda (gene 1): 1.2387
#> epsilon (gene 1): -0.528371
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 668.054
#> delta (gene 1): 0.861791
#> phi (cell 1): 0.856606
#> s (cell 1): 0.295293
#> nu (cell 1): 0.446113
#> theta (batch 1): 0.162491
#> betas:    0.2585   0.0202   0.9159   1.1612   0.4752   0.9552  -0.1269   0.3601   0.0840  -0.1111  -0.6265  -0.1834
#> 
#> sigma: 0.146253
#> lambda (gene 1): 1.79816
#> epsilon (gene 1): -0.0289238
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 701.03
#> delta (gene 1): 0.546107
#> phi (cell 1): 0.89453
#> s (cell 1): 0.461921
#> nu (cell 1): 0.421711
#> theta (batch 1): 0.248041
#> betas:    0.7327   0.0076   0.4743   0.5720   0.7104   0.8411  -0.0623  -0.0009   0.0285  -0.3915  -0.7019  -0.5368
#> 
#> sigma: 0.218536
#> lambda (gene 1): 0.574423
#> epsilon (gene 1): -0.325291
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 605.917
#> delta (gene 1): 0.659698
#> phi (cell 1): 0.952779
#> s (cell 1): 0.417853
#> nu (cell 1): 0.419417
#> theta (batch 1): 0.0452632
#> betas:    0.7454  -0.0616   0.2523   1.1284   0.7075   0.4611   0.2916   0.1131   0.1178  -0.4788  -0.1007  -0.5265
#> 
#> sigma: 0.179051
#> lambda (gene 1): 0.526523
#> epsilon (gene 1): -0.156585
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 751.836
#> delta (gene 1): 0.658056
#> phi (cell 1): 0.912268
#> s (cell 1): 0.779427
#> nu (cell 1): 0.43532
#> theta (batch 1): 0.0518544
#> betas:    0.1100   0.0575   0.8202   1.0275   0.9717   0.3662   0.3506   0.3074  -0.3025   0.2643  -0.7614  -0.5705
#> 
#> sigma: 0.167281
#> lambda (gene 1): 0.555211
#> epsilon (gene 1): -0.308636
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 617.55
#> delta (gene 1): 0.798606
#> phi (cell 1): 0.772254
#> s (cell 1): 0.417955
#> nu (cell 1): 0.428225
#> theta (batch 1): 0.119685
#> betas:    1.2302  -0.1146   1.2410   0.5183  -0.2524   1.1016  -0.1050  -0.0235   0.3523  -0.3388  -0.4532  -0.4039
#> 
#> sigma: 0.183694
#> lambda (gene 1): 0.926375
#> epsilon (gene 1): -0.115
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 638.342
#> delta (gene 1): 0.542656
#> phi (cell 1): 0.85102
#> s (cell 1): 0.488852
#> nu (cell 1): 0.408056
#> theta (batch 1): 0.174567
#> betas:    0.3643  -0.0844   0.5645   0.6456   1.1666   0.3093   0.6225   0.4334  -0.5025   0.8128  -0.6276  -0.0450
#> 
#> sigma: 0.170184
#> lambda (gene 1): 0.424617
#> epsilon (gene 1): -0.519191
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 807.848
#> delta (gene 1): 0.730529
#> phi (cell 1): 0.77129
#> s (cell 1): 0.502882
#> nu (cell 1): 0.417021
#> theta (batch 1): 0.0113039
#> betas:    0.7984  -0.1058   0.0383   0.7754   1.0679   0.5254   0.2750   0.3934  -0.2221  -0.0268  -0.1975  -0.2233
#> 
#> sigma: 0.124474
#> lambda (gene 1): 2.03263
#> epsilon (gene 1): -0.115709
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 584.137
#> delta (gene 1): 0.648045
#> phi (cell 1): 0.735078
#> s (cell 1): 0.545499
#> nu (cell 1): 0.427705
#> theta (batch 1): 0.0348765
#> betas:    0.6679  -0.0620   0.9496   0.5583   0.3567   0.9294   0.5734  -0.4897   0.5535  -0.4455  -0.3636  -0.3866
#> 
#> sigma: 0.207139
#> lambda (gene 1): 0.717437
#> epsilon (gene 1): -0.133771
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 744.34
#> delta (gene 1): 0.803294
#> phi (cell 1): 0.689415
#> s (cell 1): 0.308132
#> nu (cell 1): 0.429721
#> theta (batch 1): 0.0462117
#> betas:    1.0947  -0.2438   0.1254   1.1949   0.2028   0.7361   0.3701   0.5150  -0.0258  -0.1347   0.5428  -0.3403
#> 
#> sigma: 0.170655
#> lambda (gene 1): 2.9401
#> epsilon (gene 1): 0.0209327
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 657.994
#> delta (gene 1): 0.541073
#> phi (cell 1): 0.709293
#> s (cell 1): 0.433042
#> nu (cell 1): 0.43609
#> theta (batch 1): 0.0148199
#> betas:    0.6765  -0.0370   0.3915   0.8913   0.8406   0.2506   0.5428   0.2233  -0.3686  -0.2153  -0.2763  -0.5984
#> 
#> sigma: 0.109391
#> lambda (gene 1): 2.03856
#> epsilon (gene 1): -0.272638
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 618.317
#> delta (gene 1): 0.713492
#> phi (cell 1): 0.787768
#> s (cell 1): 0.457754
#> nu (cell 1): 0.422392
#> theta (batch 1): 0.0265704
#> betas:    0.1549   0.0023   1.1077   0.7475   0.9836   0.9638   0.0062   0.5182   0.0876  -0.3821   0.3419  -0.8927
#> 
#> sigma: 0.219962
#> lambda (gene 1): 0.346489
#> epsilon (gene 1): -0.297186
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 596.539
#> delta (gene 1): 0.586736
#> phi (cell 1): 0.845719
#> s (cell 1): 0.516337
#> nu (cell 1): 0.42848
#> theta (batch 1): 0.0254691
#> betas:    0.7315  -0.0572   0.7231   1.1478   0.2010   0.5731   0.1149   0.6805  -0.6112   0.0804  -0.2470  -0.3940
#> 
#> sigma: 0.15728
#> lambda (gene 1): 1.14749
#> epsilon (gene 1): -0.413778
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 815.454
#> delta (gene 1): 0.781277
#> phi (cell 1): 0.769123
#> s (cell 1): 0.353913
#> nu (cell 1): 0.42446
#> theta (batch 1): 0.0655298
#> betas:    0.8388  -0.1291   0.6481   1.2531  -0.0336   0.9163   0.2146   0.5844  -0.1749  -0.3697   0.5229  -0.8778
#> 
#> sigma: 0.203766
#> lambda (gene 1): 3.72956
#> epsilon (gene 1): -0.0213561
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 537.407
#> delta (gene 1): 1.0442
#> phi (cell 1): 0.801202
#> s (cell 1): 0.750915
#> nu (cell 1): 0.419262
#> theta (batch 1): 0.0980035
#> betas:    0.3993   0.0206   0.2208   1.2743   0.8134   0.5710   0.0494   0.5206  -0.3697   0.0034  -0.6189  -0.5020
#> 
#> sigma: 0.183313
#> lambda (gene 1): 2.33318
#> epsilon (gene 1): 0.248775
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 852.336
#> delta (gene 1): 0.484171
#> phi (cell 1): 0.729455
#> s (cell 1): 0.6659
#> nu (cell 1): 0.420051
#> theta (batch 1): 0.102491
#> betas:    0.9420  -0.0470   0.2040   1.4136   0.2569   0.4297   0.4469  -0.3106   0.4069  -0.6358  -0.1515  -0.3751
#> 
#> sigma: 0.146706
#> lambda (gene 1): 2.70615
#> epsilon (gene 1): -0.711917
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 611.624
#> delta (gene 1): 0.877939
#> phi (cell 1): 0.721129
#> s (cell 1): 0.204275
#> nu (cell 1): 0.433217
#> theta (batch 1): 0.656369
#> betas:    0.7313  -0.0437   0.6352   0.4567   0.8925   0.4312   0.1607   0.3855  -0.4191   0.0355  -0.0610  -0.8750
#> 
#> sigma: 0.15395
#> lambda (gene 1): 0.835565
#> epsilon (gene 1): -0.120502
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 531.007
#> delta (gene 1): 0.639782
#> phi (cell 1): 0.67829
#> s (cell 1): 0.356782
#> nu (cell 1): 0.433061
#> theta (batch 1): 0.472847
#> betas:    1.1202  -0.0895   0.3863   0.4461   0.9579   0.1133   0.2195   0.2216  -0.2576  -0.1095  -0.5912  -0.5062
#> 
#> sigma: 0.126169
#> lambda (gene 1): 1.21186
#> epsilon (gene 1): -0.151013
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 558.699
#> delta (gene 1): 0.794948
#> phi (cell 1): 0.753219
#> s (cell 1): 0.265927
#> nu (cell 1): 0.422647
#> theta (batch 1): 0.177234
#> betas:    1.0518  -0.0566   0.4431   0.6008   0.6537   0.5838   0.1381   0.2084  -0.2562  -0.4184   0.0348  -1.0744
#> 
#> sigma: 0.165771
#> lambda (gene 1): 0.972737
#> epsilon (gene 1): -0.0855882
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 798.244
#> delta (gene 1): 0.652164
#> phi (cell 1): 0.610732
#> s (cell 1): 0.318714
#> nu (cell 1): 0.443054
#> theta (batch 1): 0.379755
#> betas:    1.1597  -0.1288   0.4685   0.4554   0.5631   0.3442   0.7245  -0.2744  -0.0161   0.1622  -0.6779  -0.1337
#> 
#> sigma: 0.217646
#> lambda (gene 1): 0.924938
#> epsilon (gene 1): -0.293414
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 770.927
#> delta (gene 1): 0.606397
#> phi (cell 1): 0.666933
#> s (cell 1): 0.615234
#> nu (cell 1): 0.419448
#> theta (batch 1): 0.180848
#> betas:    1.5193  -0.2263   0.0846   0.5960   0.2943   0.8835   0.2141  -0.2471   0.5634  -0.6537   0.0432  -0.1283
#> 
#> sigma: 0.160921
#> lambda (gene 1): 0.864572
#> epsilon (gene 1): -0.181655
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 555.468
#> delta (gene 1): 0.739169
#> phi (cell 1): 0.67927
#> s (cell 1): 0.465385
#> nu (cell 1): 0.43636
#> theta (batch 1): 0.135353
#> betas:    1.5967  -0.1870   0.5251   0.7718   0.2080   0.2359   0.1731   0.1865  -0.1823  -0.4768  -0.0881   0.2632
#> 
#> sigma: 0.148828
#> lambda (gene 1): 0.300736
#> epsilon (gene 1): -0.117518
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 622.439
#> delta (gene 1): 0.66042
#> phi (cell 1): 0.769584
#> s (cell 1): 0.365902
#> nu (cell 1): 0.418588
#> theta (batch 1): 0.108811
#> betas:    1.0812  -0.0545   0.8414   0.1820   0.9594   0.2349   0.0870  -0.0412   0.0008  -0.4371  -0.4502  -0.6185
#> 
#> sigma: 0.148295
#> lambda (gene 1): 0.561665
#> epsilon (gene 1): -0.170086
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 787.936
#> delta (gene 1): 0.794249
#> phi (cell 1): 0.819961
#> s (cell 1): 0.261442
#> nu (cell 1): 0.424894
#> theta (batch 1): 0.62181
#> betas:    0.4020  -0.0894   0.6064   1.2660   0.6633   0.5466   0.4939   0.5647  -0.1693  -0.0185  -0.0266  -0.1053
#> 
#> sigma: 0.180997
#> lambda (gene 1): 1.33788
#> epsilon (gene 1): 0.00383399
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 695.772
#> delta (gene 1): 0.733895
#> phi (cell 1): 0.782725
#> s (cell 1): 0.371522
#> nu (cell 1): 0.416376
#> theta (batch 1): 0.466274
#> betas:    1.0727  -0.1595   1.0302   0.7699   0.3457   0.4499   0.5537  -0.0797   0.1948  -0.0807  -0.2532  -0.0128
#> 
#> sigma: 0.15911
#> lambda (gene 1): 1.26803
#> epsilon (gene 1): -0.173402
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 639.591
#> delta (gene 1): 0.540751
#> phi (cell 1): 0.923385
#> s (cell 1): 0.319557
#> nu (cell 1): 0.42404
#> theta (batch 1): 0.305482
#> betas:    0.8226  -0.0564   0.5343   0.5127   1.3460   0.0906   0.4525   0.2802  -0.1936  -0.1312  -0.4513  -0.4140
#> 
#> sigma: 0.142475
#> lambda (gene 1): 1.07163
#> epsilon (gene 1): -0.416119
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 821.651
#> delta (gene 1): 0.574637
#> phi (cell 1): 0.982537
#> s (cell 1): 0.349542
#> nu (cell 1): 0.414259
#> theta (batch 1): 0.0953367
#> betas:    0.7656  -0.0021   0.2683   1.3498   0.5228   0.2572   0.4920  -0.6424   0.7025  -0.6206  -0.6832  -0.3895
#> 
#> sigma: 0.173652
#> lambda (gene 1): 1.03952
#> epsilon (gene 1): -0.390839
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 860.928
#> delta (gene 1): 0.535335
#> phi (cell 1): 0.912639
#> s (cell 1): 0.374526
#> nu (cell 1): 0.433273
#> theta (batch 1): 0.109885
#> betas:    1.1239  -0.1967   0.5709   0.5124   0.8356   0.1087   0.5979   0.2433   0.0997  -0.0750  -0.2864   0.2113
#> 
#> sigma: 0.171723
#> lambda (gene 1): 0.794332
#> epsilon (gene 1): -0.345231
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 686.261
#> delta (gene 1): 0.726599
#> phi (cell 1): 0.739586
#> s (cell 1): 0.645649
#> nu (cell 1): 0.428182
#> theta (batch 1): 0.517727
#> betas:    1.0873  -0.1398   0.5317   0.3070   0.9414   0.0495   0.6702   0.1670  -0.2333  -0.2382   0.0415  -0.1658
#> 
#> sigma: 0.134211
#> lambda (gene 1): 1.98777
#> epsilon (gene 1): -0.129784
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 18980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 567.223
#> delta (gene 1): 0.543203
#> phi (cell 1): 0.727776
#> s (cell 1): 0.624792
#> nu (cell 1): 0.438622
#> theta (batch 1): 0.262098
#> betas:    0.8129  -0.1007   0.6203   0.7774   0.5053   0.6573   0.0835   0.4684  -0.1182  -0.3348  -0.0776  -0.2183
#> 
#> sigma: 0.202154
#> lambda (gene 1): 2.34885
#> epsilon (gene 1): -0.324062
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19000 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 673.083
#> delta (gene 1): 0.821707
#> phi (cell 1): 0.778205
#> s (cell 1): 0.290093
#> nu (cell 1): 0.435395
#> theta (batch 1): 0.241609
#> betas:    0.9169  -0.0902   0.3708   0.7054   0.6400   0.1662   0.6397   0.0576  -0.2609   0.2179  -0.7137  -0.0621
#> 
#> sigma: 0.154558
#> lambda (gene 1): 1.05246
#> epsilon (gene 1): -0.0330112
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19020 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 481.345
#> delta (gene 1): 0.909663
#> phi (cell 1): 0.614952
#> s (cell 1): 0.576856
#> nu (cell 1): 0.417565
#> theta (batch 1): 0.207082
#> betas:    1.1462  -0.0959   0.6268   0.3360   0.4043   0.6933   0.1772  -0.0498  -0.1256  -0.1773  -0.2267  -0.4552
#> 
#> sigma: 0.132528
#> lambda (gene 1): 0.358903
#> epsilon (gene 1): -0.0369224
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19040 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 512.24
#> delta (gene 1): 0.913321
#> phi (cell 1): 0.588211
#> s (cell 1): 0.226776
#> nu (cell 1): 0.434508
#> theta (batch 1): 0.342081
#> betas:    1.2277  -0.0841   0.4247   0.7542   0.3812   0.3534   0.1066   0.1050   0.3214  -0.9667  -0.3044  -0.3353
#> 
#> sigma: 0.146678
#> lambda (gene 1): 1.09877
#> epsilon (gene 1): 0.297695
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19060 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 622.501
#> delta (gene 1): 0.665054
#> phi (cell 1): 0.66917
#> s (cell 1): 0.303535
#> nu (cell 1): 0.432182
#> theta (batch 1): 0.209834
#> betas:    0.3977  -0.0068   0.6413   1.0158   0.4360   1.1175   0.0672   0.0490   0.4026  -0.5737  -0.2702  -0.1699
#> 
#> sigma: 0.141983
#> lambda (gene 1): 0.363738
#> epsilon (gene 1): -0.210541
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19080 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.722
#> delta (gene 1): 0.575269
#> phi (cell 1): 0.805793
#> s (cell 1): 0.72527
#> nu (cell 1): 0.419401
#> theta (batch 1): 0.608763
#> betas:    0.3591  -0.0168   0.9469   0.8112   0.5731   0.8021   0.3615   0.0877   0.3524  -0.6173   0.0885  -0.6480
#> 
#> sigma: 0.155532
#> lambda (gene 1): 2.82355
#> epsilon (gene 1): -0.32616
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19100 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 787.518
#> delta (gene 1): 0.656049
#> phi (cell 1): 0.916149
#> s (cell 1): 0.363112
#> nu (cell 1): 0.43218
#> theta (batch 1): 0.135858
#> betas:    0.7007   0.1146   0.3624   1.1619   0.3482   0.3071   0.5245  -0.3947  -0.2335  -0.5348  -0.8321  -0.8667
#> 
#> sigma: 0.189804
#> lambda (gene 1): 0.911119
#> epsilon (gene 1): -0.249165
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19120 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 807.666
#> delta (gene 1): 0.617565
#> phi (cell 1): 0.808877
#> s (cell 1): 0.465389
#> nu (cell 1): 0.423251
#> theta (batch 1): 0.0707001
#> betas:    0.7343  -0.0557   1.2506   0.7587   0.5904   0.2261   0.5592  -0.1633   0.4079  -0.7014  -0.0507  -0.3934
#> 
#> sigma: 0.15502
#> lambda (gene 1): 0.412419
#> epsilon (gene 1): -0.244013
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19140 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 571.701
#> delta (gene 1): 0.677311
#> phi (cell 1): 0.734204
#> s (cell 1): 0.579763
#> nu (cell 1): 0.422793
#> theta (batch 1): 0.111998
#> betas:    0.6247   0.0208   0.2219   1.2058   0.5776   0.3022   0.3237  -0.2513   0.6210  -1.0619  -0.1819  -0.4330
#> 
#> sigma: 0.155983
#> lambda (gene 1): 1.07938
#> epsilon (gene 1): -0.192362
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19160 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 930.364
#> delta (gene 1): 0.782377
#> phi (cell 1): 0.754172
#> s (cell 1): 0.534508
#> nu (cell 1): 0.425778
#> theta (batch 1): 0.0306169
#> betas:    0.6196  -0.0616   0.6601   0.6981   0.5416   1.0142   0.1865  -0.0225  -0.0579   0.2095  -0.5119  -0.4264
#> 
#> sigma: 0.144146
#> lambda (gene 1): 1.41169
#> epsilon (gene 1): -0.115686
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19180 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 635.132
#> delta (gene 1): 0.660788
#> phi (cell 1): 0.855142
#> s (cell 1): 0.394441
#> nu (cell 1): 0.423573
#> theta (batch 1): 0.0739705
#> betas:    0.7964  -0.1858   1.3140   0.4913   0.7134   0.5378   0.4566   0.1490   0.4150   0.0368  -0.2693   0.5659
#> 
#> sigma: 0.177786
#> lambda (gene 1): 1.3598
#> epsilon (gene 1): -0.28874
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19200 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 672.311
#> delta (gene 1): 0.688103
#> phi (cell 1): 0.894618
#> s (cell 1): 0.507511
#> nu (cell 1): 0.439161
#> theta (batch 1): 0.0427969
#> betas:    0.9915  -0.2484   0.2321   0.8143   1.5449  -0.0641   0.6268   0.3640   0.4052   0.0061   0.2196   0.0783
#> 
#> sigma: 0.207987
#> lambda (gene 1): 1.78563
#> epsilon (gene 1): -0.269889
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19220 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 691.837
#> delta (gene 1): 0.626975
#> phi (cell 1): 0.866558
#> s (cell 1): 0.334112
#> nu (cell 1): 0.409013
#> theta (batch 1): 0.169588
#> betas:    0.1235   0.0365   0.7001   0.8960   0.4885   1.1808   0.0697  -0.0126   0.6104  -0.4528  -0.6566   0.1898
#> 
#> sigma: 0.174471
#> lambda (gene 1): 1.50602
#> epsilon (gene 1): -0.336133
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19240 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 743.524
#> delta (gene 1): 0.480615
#> phi (cell 1): 0.768495
#> s (cell 1): 0.444265
#> nu (cell 1): 0.438736
#> theta (batch 1): 0.265377
#> betas:    0.9313  -0.1789   0.2060   0.9395   0.9570   0.2671   0.5398  -0.0055   0.2154   0.1571  -0.0776  -0.1359
#> 
#> sigma: 0.137904
#> lambda (gene 1): 0.33644
#> epsilon (gene 1): -0.689349
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19260 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 803.99
#> delta (gene 1): 0.698819
#> phi (cell 1): 0.686649
#> s (cell 1): 0.672858
#> nu (cell 1): 0.410413
#> theta (batch 1): 0.277887
#> betas:    1.1431  -0.2580  -0.0088   1.1587   0.5359   0.7955   0.5825   0.0498   0.2172   0.1360   0.1041   0.0789
#> 
#> sigma: 0.186256
#> lambda (gene 1): 0.510499
#> epsilon (gene 1): -0.195912
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19280 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 741.605
#> delta (gene 1): 0.605468
#> phi (cell 1): 0.70836
#> s (cell 1): 2.05129
#> nu (cell 1): 0.42472
#> theta (batch 1): 0.51828
#> betas:    0.5033  -0.0724   1.4482   0.6663   0.8714   0.3424   0.2750   0.2068   0.2802  -0.0636  -0.4064  -0.4442
#> 
#> sigma: 0.186473
#> lambda (gene 1): 1.1108
#> epsilon (gene 1): -0.256198
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19300 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 692.952
#> delta (gene 1): 0.840256
#> phi (cell 1): 0.695439
#> s (cell 1): 0.641284
#> nu (cell 1): 0.419396
#> theta (batch 1): 0.327309
#> betas:    0.6795   0.0369   0.3139   1.0662   0.4027   0.5444   0.2105  -0.4364   0.4444  -0.7537  -0.5452  -0.4954
#> 
#> sigma: 0.245097
#> lambda (gene 1): 0.662019
#> epsilon (gene 1): 0.00633092
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19320 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 633.694
#> delta (gene 1): 0.781034
#> phi (cell 1): 0.568413
#> s (cell 1): 0.413251
#> nu (cell 1): 0.426825
#> theta (batch 1): 0.131662
#> betas:    1.0781  -0.2140   1.1984   0.6294   0.5998   0.5121   0.5276  -0.2393   0.7322   0.2504  -0.6117   0.0895
#> 
#> sigma: 0.162498
#> lambda (gene 1): 2.35377
#> epsilon (gene 1): -0.154031
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19340 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 539.861
#> delta (gene 1): 0.613007
#> phi (cell 1): 0.626452
#> s (cell 1): 0.333977
#> nu (cell 1): 0.424617
#> theta (batch 1): 0.283988
#> betas:    0.9781  -0.0935   0.8337   0.1566   0.8088   0.8469  -0.1642   0.3364  -0.3677  -0.0119  -0.4477  -0.3007
#> 
#> sigma: 0.115405
#> lambda (gene 1): 0.954978
#> epsilon (gene 1): -0.233434
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19360 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 903.191
#> delta (gene 1): 0.58409
#> phi (cell 1): 0.791343
#> s (cell 1): 0.369867
#> nu (cell 1): 0.42511
#> theta (batch 1): 0.191151
#> betas:    1.4733  -0.2782   0.3454   0.5808   0.3166   0.7514   0.7090  -0.3763   0.5994  -0.1384   0.1255   0.0151
#> 
#> sigma: 0.184375
#> lambda (gene 1): 0.87927
#> epsilon (gene 1): -0.458032
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19380 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 745.117
#> delta (gene 1): 0.554938
#> phi (cell 1): 0.772923
#> s (cell 1): 0.413447
#> nu (cell 1): 0.416132
#> theta (batch 1): 0.279179
#> betas:    0.1411   0.0819   0.1952   1.4133   0.7202   0.2297   0.5988  -0.1643  -0.1982   0.2646  -0.9665  -0.3900
#> 
#> sigma: 0.173524
#> lambda (gene 1): 0.228241
#> epsilon (gene 1): -0.545696
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19400 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 801.493
#> delta (gene 1): 0.755629
#> phi (cell 1): 0.736662
#> s (cell 1): 0.572089
#> nu (cell 1): 0.424833
#> theta (batch 1): 0.233299
#> betas:    0.8284  -0.1562   0.4010   1.0650   0.2672   0.7681   0.3809   0.4528   0.0086   0.0951  -0.3273   0.2339
#> 
#> sigma: 0.1432
#> lambda (gene 1): 0.303713
#> epsilon (gene 1): -0.0987133
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19420 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 730.895
#> delta (gene 1): 0.640599
#> phi (cell 1): 0.699637
#> s (cell 1): 0.684942
#> nu (cell 1): 0.420924
#> theta (batch 1): 0.44234
#> betas:    1.0351  -0.1448   0.7830   1.1486   0.4016   0.6367   0.1725   0.4046   0.0681  -0.2410  -0.1315  -0.3432
#> 
#> sigma: 0.151087
#> lambda (gene 1): 1.14338
#> epsilon (gene 1): -0.213953
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19440 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 803.283
#> delta (gene 1): 0.604117
#> phi (cell 1): 0.798044
#> s (cell 1): 0.677116
#> nu (cell 1): 0.428196
#> theta (batch 1): 0.194048
#> betas:    0.5385  -0.0403   0.7964   0.8332   0.7913   0.2855   0.8150   0.1153  -0.0929  -0.2019  -0.3725  -0.3317
#> 
#> sigma: 0.170615
#> lambda (gene 1): 0.585059
#> epsilon (gene 1): -0.185108
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19460 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 613.958
#> delta (gene 1): 0.516183
#> phi (cell 1): 0.629469
#> s (cell 1): 0.435786
#> nu (cell 1): 0.410852
#> theta (batch 1): 0.0892783
#> betas:    0.8953  -0.0528   0.5643   0.6881   0.9107   0.3982   0.0990  -0.1369  -0.0086  -0.1203  -0.7196  -0.3456
#> 
#> sigma: 0.159162
#> lambda (gene 1): 0.803903
#> epsilon (gene 1): -0.405732
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19480 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 629.194
#> delta (gene 1): 0.724406
#> phi (cell 1): 0.592524
#> s (cell 1): 0.556507
#> nu (cell 1): 0.430651
#> theta (batch 1): 0.14594
#> betas:    0.9935  -0.1588   0.6414   0.8032   0.6698   0.1987   0.4529   0.5541  -0.4662  -0.0620   0.1320  -0.2383
#> 
#> sigma: 0.164259
#> lambda (gene 1): 0.491021
#> epsilon (gene 1): -0.0593788
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19500 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 572.544
#> delta (gene 1): 0.7362
#> phi (cell 1): 0.695175
#> s (cell 1): 0.179092
#> nu (cell 1): 0.412762
#> theta (batch 1): 0.532302
#> betas:    0.9638  -0.1282   0.8980   0.0736   1.0708   0.3361   0.4035   0.1628   0.0578  -0.1969  -0.2474  -0.2088
#> 
#> sigma: 0.143311
#> lambda (gene 1): 1.77207
#> epsilon (gene 1): -0.0583983
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19520 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 803.345
#> delta (gene 1): 0.596638
#> phi (cell 1): 0.694434
#> s (cell 1): 0.804146
#> nu (cell 1): 0.435527
#> theta (batch 1): 0.307397
#> betas:    0.3583   0.0678   0.5987   0.9948   0.5702   0.3804   0.3465  -0.0327   0.0344  -0.4615  -0.7060  -0.3710
#> 
#> sigma: 0.161946
#> lambda (gene 1): 0.232842
#> epsilon (gene 1): -0.248047
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19540 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 795.634
#> delta (gene 1): 0.705673
#> phi (cell 1): 0.773263
#> s (cell 1): 0.519977
#> nu (cell 1): 0.40889
#> theta (batch 1): 0.119382
#> betas:    0.8955  -0.1298   0.7546   0.6846   0.9627   0.0988   0.7004   0.3337  -0.3439  -0.0906   0.0343  -0.4102
#> 
#> sigma: 0.121132
#> lambda (gene 1): 0.597341
#> epsilon (gene 1): -0.0684895
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19560 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 715.072
#> delta (gene 1): 0.801972
#> phi (cell 1): 0.728526
#> s (cell 1): 0.535302
#> nu (cell 1): 0.434867
#> theta (batch 1): 0.147356
#> betas:    1.2266  -0.1296  -0.0470   0.8565   1.3899  -0.2412   0.2445   0.0958   0.1688  -0.4860  -0.3619  -0.1868
#> 
#> sigma: 0.136696
#> lambda (gene 1): 1.46477
#> epsilon (gene 1): 0.0918844
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19580 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 874.496
#> delta (gene 1): 0.887887
#> phi (cell 1): 0.753962
#> s (cell 1): 3.61154
#> nu (cell 1): 0.423039
#> theta (batch 1): 0.4723
#> betas:    1.3408  -0.0156   0.3473   0.9439   0.2714   0.2463   0.2119  -0.3875  -0.2779  -0.2625  -1.0611  -0.5706
#> 
#> sigma: 0.162238
#> lambda (gene 1): 1.05131
#> epsilon (gene 1): 0.13543
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19600 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 729.522
#> delta (gene 1): 0.737272
#> phi (cell 1): 0.643432
#> s (cell 1): 0.927657
#> nu (cell 1): 0.434015
#> theta (batch 1): 0.287603
#> betas:    0.7598  -0.1790   1.3067   0.5331   0.9934   0.5994   0.1569   0.9790  -0.3024  -0.3530   0.5049  -0.1962
#> 
#> sigma: 0.172271
#> lambda (gene 1): 0.974119
#> epsilon (gene 1): 0.132624
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19620 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 685.035
#> delta (gene 1): 0.509884
#> phi (cell 1): 0.73659
#> s (cell 1): 0.271633
#> nu (cell 1): 0.406806
#> theta (batch 1): 0.40903
#> betas:    0.6294   0.0484   0.7432   0.3456   0.2648   1.3220  -0.0564  -0.4091   0.3789  -0.6310  -1.0665   0.0238
#> 
#> sigma: 0.179954
#> lambda (gene 1): 0.517038
#> epsilon (gene 1): -0.351749
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19640 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 622.228
#> delta (gene 1): 0.705956
#> phi (cell 1): 0.707187
#> s (cell 1): 0.574874
#> nu (cell 1): 0.427387
#> theta (batch 1): 0.310361
#> betas:    0.9732  -0.1621   0.6517   0.4168   0.7832   0.4205   0.2732   0.0790   0.4705  -0.1857  -0.0940  -0.1913
#> 
#> sigma: 0.142157
#> lambda (gene 1): 1.11246
#> epsilon (gene 1): -0.262872
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19660 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 617.612
#> delta (gene 1): 0.66024
#> phi (cell 1): 0.69517
#> s (cell 1): 1.05269
#> nu (cell 1): 0.422844
#> theta (batch 1): 0.177309
#> betas:    0.5857  -0.0687   0.0996   1.1896   0.7870   0.4714   0.3046   0.1722  -0.0053  -0.1647  -0.2461  -0.5046
#> 
#> sigma: 0.136892
#> lambda (gene 1): 0.938399
#> epsilon (gene 1): -0.0894106
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19680 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 693.367
#> delta (gene 1): 0.593126
#> phi (cell 1): 0.708224
#> s (cell 1): 0.377132
#> nu (cell 1): 0.422033
#> theta (batch 1): 0.0371008
#> betas:    0.8713  -0.0283   0.7347   0.5954   0.5010   0.6583   0.0279   0.1836  -0.2930  -0.3492  -0.4010  -0.5914
#> 
#> sigma: 0.175112
#> lambda (gene 1): 0.309499
#> epsilon (gene 1): -0.230852
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19700 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 832.792
#> delta (gene 1): 0.571163
#> phi (cell 1): 0.628425
#> s (cell 1): 0.599006
#> nu (cell 1): 0.447172
#> theta (batch 1): 0.111861
#> betas:    0.9829  -0.1176   0.7553   0.6093   0.5020   0.5645   0.2361   0.4143  -0.2805   0.1126  -0.5354  -0.0496
#> 
#> sigma: 0.155819
#> lambda (gene 1): 1.34342
#> epsilon (gene 1): -0.365326
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19720 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 646.132
#> delta (gene 1): 0.641698
#> phi (cell 1): 0.67076
#> s (cell 1): 0.473449
#> nu (cell 1): 0.42895
#> theta (batch 1): 0.110654
#> betas:    1.3534  -0.1542   0.6866   0.6445   0.2673   0.4938   0.0633   0.0744  -0.0911  -0.1425  -0.2679  -0.4499
#> 
#> sigma: 0.136954
#> lambda (gene 1): 1.06529
#> epsilon (gene 1): -0.253433
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19740 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 670.721
#> delta (gene 1): 0.696297
#> phi (cell 1): 0.689653
#> s (cell 1): 0.354778
#> nu (cell 1): 0.429366
#> theta (batch 1): 0.123812
#> betas:    1.5676  -0.1693   0.7260   1.0213   0.0537   0.4970   0.2596   0.2988  -0.6587  -0.1025  -0.1276  -0.4507
#> 
#> sigma: 0.150573
#> lambda (gene 1): 0.677623
#> epsilon (gene 1): -0.11803
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19760 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 712.088
#> delta (gene 1): 0.629502
#> phi (cell 1): 0.711922
#> s (cell 1): 0.61965
#> nu (cell 1): 0.411538
#> theta (batch 1): 0.373213
#> betas:    0.3715  -0.0389   0.6022   0.8545   0.9171   0.6563   0.0595   0.2008   0.4418  -0.5151  -0.3144   0.0793
#> 
#> sigma: 0.112214
#> lambda (gene 1): 0.4151
#> epsilon (gene 1): -0.178771
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19780 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 601.828
#> delta (gene 1): 0.670209
#> phi (cell 1): 0.731268
#> s (cell 1): 0.484646
#> nu (cell 1): 0.404841
#> theta (batch 1): 0.165268
#> betas:    0.8671  -0.0902   0.4138   0.6631   0.8919   0.5623  -0.1969   0.6043  -0.0338  -0.1224  -0.6445   0.0236
#> 
#> sigma: 0.194351
#> lambda (gene 1): 0.744693
#> epsilon (gene 1): -0.160934
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19800 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 807.186
#> delta (gene 1): 0.945372
#> phi (cell 1): 0.780236
#> s (cell 1): 1.19158
#> nu (cell 1): 0.418253
#> theta (batch 1): 0.382286
#> betas:    1.0819  -0.1340   0.1171   0.4850   1.1836   0.0913   0.3715   0.2616  -0.1608  -0.1675  -0.1293  -0.2908
#> 
#> sigma: 0.103047
#> lambda (gene 1): 0.670983
#> epsilon (gene 1): 0.139854
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19820 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 844.493
#> delta (gene 1): 0.82057
#> phi (cell 1): 0.694072
#> s (cell 1): 0.35728
#> nu (cell 1): 0.419217
#> theta (batch 1): 0.598372
#> betas:    0.3507  -0.0421   0.2747   0.9358   1.4136   0.4968   0.0376   0.6049   0.0845  -0.2647  -0.1580  -0.2026
#> 
#> sigma: 0.15936
#> lambda (gene 1): 0.702407
#> epsilon (gene 1): -0.0107474
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19840 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 537.789
#> delta (gene 1): 0.889667
#> phi (cell 1): 0.689829
#> s (cell 1): 0.278642
#> nu (cell 1): 0.423322
#> theta (batch 1): 0.109554
#> betas:    0.8382  -0.1699   1.0302   0.6885   0.8026   0.8495   0.2257   0.2686   0.1361  -0.1222   0.1273  -0.0758
#> 
#> sigma: 0.146662
#> lambda (gene 1): 0.770347
#> epsilon (gene 1): 0.066926
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19860 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 696.056
#> delta (gene 1): 0.969199
#> phi (cell 1): 0.598482
#> s (cell 1): 0.399195
#> nu (cell 1): 0.409437
#> theta (batch 1): 0.179885
#> betas:    0.9115  -0.0818   0.7190   0.5524   0.7339   0.0808   0.8298  -0.0270  -0.2603   0.3023  -0.9506  -0.1971
#> 
#> sigma: 0.174772
#> lambda (gene 1): 0.82547
#> epsilon (gene 1): 0.228972
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19880 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 735.8
#> delta (gene 1): 0.434154
#> phi (cell 1): 0.622039
#> s (cell 1): 0.885941
#> nu (cell 1): 0.415538
#> theta (batch 1): 0.420152
#> betas:    1.1602  -0.1780   0.8314   0.6433   0.6443  -0.0579   0.7975   0.1403   0.0175  -0.2488   0.0073  -0.1081
#> 
#> sigma: 0.141154
#> lambda (gene 1): 0.513544
#> epsilon (gene 1): -0.604771
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19900 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 701.063
#> delta (gene 1): 0.709701
#> phi (cell 1): 0.705425
#> s (cell 1): 0.317211
#> nu (cell 1): 0.431867
#> theta (batch 1): 0.247056
#> betas:    0.5567   0.0062   0.5399   0.8166   0.7331   0.6889   0.1489   0.2656  -0.3467  -0.3623  -0.3753  -0.5710
#> 
#> sigma: 0.14769
#> lambda (gene 1): 1.52646
#> epsilon (gene 1): 0.0343678
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19920 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 891.96
#> delta (gene 1): 0.787382
#> phi (cell 1): 0.856783
#> s (cell 1): 0.43976
#> nu (cell 1): 0.430787
#> theta (batch 1): 0.17366
#> betas:    0.9390  -0.1821   0.6638   0.6278   0.2215   1.2388   0.2240   0.2651  -0.0787   0.1612  -0.3667   0.3516
#> 
#> sigma: 0.206091
#> lambda (gene 1): 1.76554
#> epsilon (gene 1): 0.0340654
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19940 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 733.874
#> delta (gene 1): 0.74082
#> phi (cell 1): 0.802076
#> s (cell 1): 0.683857
#> nu (cell 1): 0.408189
#> theta (batch 1): 0.108299
#> betas:    0.7261  -0.0676   0.8137   0.9536   0.4785   0.7251   0.1658   0.3071  -0.3950   0.1008  -0.6194  -0.2136
#> 
#> sigma: 0.127034
#> lambda (gene 1): 0.393947
#> epsilon (gene 1): 0.0416183
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19960 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 634.411
#> delta (gene 1): 0.456463
#> phi (cell 1): 0.747775
#> s (cell 1): 0.280991
#> nu (cell 1): 0.412054
#> theta (batch 1): 0.17293
#> betas:    2.1193  -0.3035  -0.1349   0.6143   0.7955   0.0122   0.1653   0.1640   0.0508  -0.2068  -0.4311   0.2118
#> 
#> sigma: 0.190779
#> lambda (gene 1): 0.887112
#> epsilon (gene 1): -0.484822
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#> -----------------------------------------------------
#> Iteration 19980 out of 20000 has been completed.
#> -----------------------------------------------------
#> Current draws for selected parameters are displayed below.
#> mu (gene 1): 932.122
#> delta (gene 1): 0.837307
#> phi (cell 1): 0.812771
#> s (cell 1): 0.394414
#> nu (cell 1): 0.416913
#> theta (batch 1): 0.42181
#> betas:    0.3850   0.0064   0.7165   0.7733   0.2240   0.8927   0.3893  -0.0333  -0.0230  -0.1289  -0.3126  -0.5920
#> 
#> sigma: 0.110629
#> lambda (gene 1): 0.810266
#> epsilon (gene 1): -0.0595287
#> -----------------------------------------------------
#> Current proposal variances for Metropolis Hastings updates (log-scale).
#> LSmu (gene 1): -2.4
#> LSdelta (gene 1): -1.6
#> LSphi: 7.2
#> LSnu (cell 1): -5.51822
#> LStheta (batch 1): -1
#>  
#> -----------------------------------------------------
#> -----------------------------------------------------
#> All 20000 MCMC iterations have been completed.
#> -----------------------------------------------------
#> -----------------------------------------------------
#>  
#> -----------------------------------------------------
#> Please see below a summary of the overall acceptance rates.
#> -----------------------------------------------------
#>  
#> Minimum acceptance rate among mu[i]'s: 0.206267
#> Average acceptance rate among mu[i]'s: 0.433219
#> Maximum acceptance rate among mu[i]'s: 0.5622
#>  
#>  
#> Minimum acceptance rate among delta[i]'s: 0.3644
#> Average acceptance rate among delta[i]'s: 0.44467
#> Maximum acceptance rate among delta[i]'s: 0.5166
#>  
#>  
#> Acceptance rate for phi (joint): 0.497533
#>  
#>  
#> Minimum acceptance rate among nu[j]'s: 0.3722
#> Average acceptance rate among nu[j]'s: 0.441709
#> Maximum acceptance rate among nu[j]'s: 0.509133
#>  
#>  
#> Minimum acceptance rate among theta[k]'s: 0.444
#> Average acceptance rate among theta[k]'s: 0.444
#> Maximum acceptance rate among theta[k]'s: 0.444
#>  
#> -----------------------------------------------------
#>  
#> -------------------------------------------------------------
#> MCMC running time 
#> -------------------------------------------------------------
#> user: 5.236
#> system: 0.177
#> elapsed: 5.414
#> -------------------------------------------------------------
#> Output 
#> -------------------------------------------------------------
#> -------------------------------------------------------------
#> BASiCS version 2.25.0 : 
#> vertical integration (spikes case) 
#> -------------------------------------------------------------
params
#> A Params object of class BASiCSParams 
#> Parameters can be (estimable) or [not estimable], 'Default' or  'NOT DEFAULT' 
#> Secondary parameters are usually set during simulation
#> 
#> Global: 
#> (GENES)  (CELLS)   [SEED] 
#>      90       30   447753 
#> 
#> 7 additional parameters 
#> 
#> Genes: 
#> 
#> (PARAMS)
#> data.frame (90 x 2) with columns: Mean, Delta 
#>                Mean     Delta
#> Gene_0011 695.35308 0.6900794
#> Gene_0012  12.29801 9.0372164
#> Gene_0013  98.06908 2.9618157
#> Gene_0014  91.65091 5.5753079
#> # ... with 86 more rows
#> 
#> Cells: 
#> 
#> (PARAMS)
#> data.frame (30 x 2) with columns: Phi, S 
#>                       Phi         S
#> Cell_001_Batch1 0.7575119 0.4701208
#> Cell_002_Batch1 1.7109542 0.2195495
#> Cell_003_Batch1 0.8980183 0.4876586
#> Cell_004_Batch1 0.5180598 0.7561904
#> # ... with 26 more rows
#> 
#> Batches: 
#>     (BATCHES)  (BATCH CELLS) 
#>             1             30 
#> 
#> Spike-ins: 
#>         (NUMBER)           (MEANS) 
#>               10  700.329643969687 
#> 
#> Variability: 
#>           (THETA) 
#> 0.196620256627407 
#> 
# }
```
