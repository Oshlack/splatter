# The LunParams class

S4 class that holds parameters for the Lun simulation.

## Parameters

The Lun simulation uses the following parameters:

- `nGenes`:

  The number of genes to simulate.

- `nCells`:

  The number of cells to simulate.

- `[nGroups]`:

  The number of groups to simulate.

- `[groupCells]`:

  Vector giving the number of cells in each simulation group/path.

- `[seed]`:

  Seed to use for generating random numbers.

- *Mean parameters*:

  `[mean.shape]`

  : Shape parameter for the mean gamma distribution.

  `[mean.rate]`

  : Rate parameter for the mean gamma distribution.

- *Counts parameters*:

  `[count.disp]`

  : The dispersion parameter for the counts negative binomial
    distribution.

- *Differential expression parameters*:

  `[de.nGenes]`

  : The number of genes that are differentially expressed in each group

  `[de.upProp]`

  : The proportion of differentially expressed genes that are
    up-regulated in each group

  `[de.upFC]`

  : The fold change for up-regulated genes

  `[de.downFC]`

  : The fold change for down-regulated genes

The parameters not shown in brackets can be estimated from real data
using
[`lunEstimate`](http://oshlacklab.com/splatter/reference/lunEstimate.md).
For details of the Lun simulation see
[`lunSimulate`](http://oshlacklab.com/splatter/reference/lunSimulate.md).
