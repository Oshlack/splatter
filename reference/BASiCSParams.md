# The BASiCSParams class

S4 class that holds parameters for the BASiCS simulation.

## Parameters

The BASiCS simulation uses the following parameters:

- `nGenes`:

  The number of genes to simulate.

- `nCells`:

  The number of cells to simulate.

- `[seed]`:

  Seed to use for generating random numbers.

- *Batch parameters*:

  `nBatches`

  : Number of batches to simulate.

  `batchCells`

  : Number of cells in each batch.

- *Gene parameters*:

  `gene.params`

  : A `data.frame` containing gene parameters with two columns: `Mean`
    (mean expression for each biological gene) and `Delta` (cell-to-cell
    heterogeneity for each biological gene).

- *Spike-in parameters*:

  `nSpikes`

  : The number of spike-ins to simulate.

  `spike.means`

  : Input molecules for each spike-in.

- *Cell parameters*:

  `cell.params`

  : A `data.frame` containing gene parameters with two columns: `Phi`
    (mRNA content factor for each cell, scaled to sum to the number of
    cells in each batch) and `S` (capture efficient for each cell).

- *Variability parameters*:

  `theta`

  : Technical variability parameter for each batch.

The parameters not shown in brackets can be estimated from real data
using
[`BASiCSEstimate`](http://oshlacklab.com/splatter/reference/BASiCSEstimate.md).
For details of the BASiCS simulation see
[`BASiCSSimulate`](http://oshlacklab.com/splatter/reference/BASiCSSimulate.md).
