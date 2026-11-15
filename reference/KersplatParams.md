# The KersplatParams class

S4 class that holds parameters for the Kersplat simulation.

## Parameters

The Kersplat simulation uses the following parameters:

- `nGenes`:

  The number of genes to simulate.

- `nCells`:

  The number of cells to simulate.

- `[seed]`:

  Seed to use for generating random numbers.

- *Mean parameters*:

  `mean.shape`

  : Shape parameter for the mean gamma distribution.

  `mean.rate`

  : Rate parameter for the mean gamma distribution.

  `mean.outProb`

  : Probability that a gene is an expression outlier.

  `mean.outFacLoc`

  : Location (meanlog) parameter for the expression outlier factor
    log-normal distribution.

  `mean.outFacScale`

  : Scale (sdlog) parameter for the expression outlier factor log-normal
    distribution.

  `mean.dens`

  : [`density`](https://rdrr.io/r/stats/density.html) object describing
    the log gene mean density.

  `[mean.method]`

  : Method to use for simulating gene means. Either "fit" to sample from
    a gamma distribution (with expression outliers) or "density" to
    sample from the provided density object.

  `[mean.values]`

  : Vector of means for each gene.

- *Biological Coefficient of Variation parameters*:

  `bcv.common`

  : Underlying common dispersion across all genes.

  `[bcv.df]`

  : Degrees of Freedom for the BCV inverse chi-squared distribution.

- *Network parameters*:

  `[network.graph]`

  : Graph containing the gene network.

  `[network.nRegs]`

  : Number of regulators in the network.

- *Paths parameters*:

  `[paths.programs]`

  : Number of expression programs.

  `[paths.design]`

  : data.frame describing path structure. See
    [`kersplatSimPaths`](http://oshlacklab.com/splatter/reference/kersplatSimPaths.md)
    for details.

- *Library size parameters*:

  `lib.loc`

  : Location (meanlog) parameter for the library size log-normal
    distribution, or mean parameter if a normal distribution is used.

  `lib.scale`

  : Scale (sdlog) parameter for the library size log-normal
    distribution, or sd parameter if a normal distribution is used.

  `lib.dens`

  : [`density`](https://rdrr.io/r/stats/density.html) object describing
    the library size density.

  `[lib.method]`

  : Method to use for simulating library sizes. Either "fit" to sample
    from a log-normal distribution or "density" to sample from the
    provided density object.

- *Design parameters*:

  `[cells.design]`

  : data.frame describing cell structure. See
    [`kersplatSimCellMeans`](http://oshlacklab.com/splatter/reference/kersplatSimCellMeans.md)
    for details.

- *Doublet parameters*:

  `[doublet.prop]`

  : Proportion of cells that are doublets.

- *Ambient parameters*:

  `[ambient.scale]`

  : Scaling factor for the library size log-normal distribution when
    generating ambient library sizes.

  `[ambient.nEmpty]`

  : Number of empty cells to simulate.

The parameters not shown in brackets can be estimated from real data
using
[`kersplatEstimate`](http://oshlacklab.com/splatter/reference/kersplatEstimate.md).
For details of the Kersplat simulation see
[`kersplatSimulate`](http://oshlacklab.com/splatter/reference/kersplatSimulate.md).
