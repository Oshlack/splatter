# The Lun2Params class

S4 class that holds parameters for the Lun2 simulation.

## Parameters

The Lun2 simulation uses the following parameters:

- `nGenes`:

  The number of genes to simulate.

- `nCells`:

  The number of cells to simulate.

- `[seed]`:

  Seed to use for generating random numbers.

- *Gene parameters*:

  `gene.params`

  : A `data.frame` containing gene parameters with two columns: `Mean`
    (mean expression for each gene) and `Disp` (dispersion for each
    gene).

  `zi.params`

  : A `data.frame` containing zero-inflated gene parameters with three
    columns: `Mean` (mean expression for each gene), `Disp` (dispersion
    for each, gene), and `Prop` (zero proportion for each gene).

- `[nPlates]`:

  The number of plates to simulate.

- *Plate parameters*:

  `plate.ingroup`

  : Character vector giving the plates considered to be part of the
    "ingroup".

  `plate.mod`

  : Plate effect modifier factor. The plate effect variance is divided
    by this value.

  `plate.var`

  : Plate effect variance.

- *Cell parameters*:

  `cell.plates`

  : Factor giving the plate that each cell comes from.

  `cell.libSizes`

  : Library size for each cell.

  `cell.libMod`

  : Modifier factor for library sizes. The library sizes are multiplied
    by this value.

- *Differential expression parameters*:

  `de.nGenes`

  : Number of differentially expressed genes.

  `de.fc`

  : Fold change for differentially expressed genes.

The parameters not shown in brackets can be estimated from real data
using
[`lun2Estimate`](http://oshlacklab.com/splatter/reference/lun2Estimate.md).
For details of the Lun2 simulation see
[`lun2Simulate`](http://oshlacklab.com/splatter/reference/lun2Simulate.md).
