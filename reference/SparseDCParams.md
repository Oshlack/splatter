# The SparseDCParams class

S4 class that holds parameters for the SparseDC simulation.

## Parameters

The SparseDC simulation uses the following parameters:

- `nGenes`:

  The number of genes to simulate in each condition.

- `nCells`:

  The number of cells to simulate.

- `[seed]`:

  Seed to use for generating random numbers.

- `markers.n`:

  Number of marker genes to simulate for each cluster.

- `markers.shared`:

  Number of marker genes for each cluster shared between conditions.
  Must be less than or equal to `markers.n`

- `[markers.same]`:

  Logical. Whether each cluster should have the same set of marker
  genes.

- `clusts.c1`:

  Numeric vector of clusters present in condition 1. The number of times
  a cluster is repeated controls the proportion of cells from that
  cluster.

- `clusts.c2`:

  Numeric vector of clusters present in condition 2. The number of times
  a cluster is repeated controls the proportion of cells from that
  cluster.

- `[mean.lower]`:

  Lower bound for cluster gene means.

- `[mean.upper]`:

  Upper bound for cluster gene means.

The parameters not shown in brackets can be estimated from real data
using
[`sparseDCEstimate`](http://oshlacklab.com/splatter/reference/sparseDCEstimate.md).
For details of the SparseDC simulation see
[`sparseDCSimulate`](http://oshlacklab.com/splatter/reference/sparseDCSimulate.md).
