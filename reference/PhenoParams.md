# The PhenoParams class

S4 class that holds parameters for the PhenoPath simulation.

## Parameters

The PhenoPath simulation uses the following parameters:

- `nGenes`:

  The number of genes to simulate.

- `nCells`:

  The number of cells to simulate.

- `[seed]`:

  Seed to use for generating random numbers.

- `[n.de]`:

  Number of genes to simulate from the differential expression regime

- `[n.pst]`:

  Number of genes to simulate from the pseudotime regime

- `[n.pst.beta]`:

  Number of genes to simulate from the pseudotime + beta interactions
  regime

- `[n.de.pst.beta]`:

  Number of genes to simulate from the differential expression +
  pseudotime + interactions regime

The parameters not shown in brackets can be estimated from real data
using
[`phenoEstimate`](http://oshlacklab.com/splatter/reference/phenoEstimate.md).
For details of the PhenoPath simulation see
[`phenoSimulate`](http://oshlacklab.com/splatter/reference/phenoSimulate.md).
