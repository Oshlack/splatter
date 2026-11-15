# The ZINBParams class

S4 class that holds parameters for the ZINB-WaVE simulation.

## Parameters

The ZINB-WaVE simulation uses the following parameters:

- `nGenes`:

  The number of genes to simulate.

- `nCells`:

  The number of cells to simulate.

- `[seed]`:

  Seed to use for generating random numbers.

- `model`:

  Object describing a ZINB model.

The majority of the parameters for this simulation are stored in a
[`ZinbModel`](https://rdrr.io/pkg/zinbwave/man/ZinbModel-class.html)
object. Please refer to the documentation for this class and its
constructor([`zinbModel`](https://rdrr.io/pkg/zinbwave/man/zinbModel.html))
for details about all the parameters.

The parameters not shown in brackets can be estimated from real data
using
[`zinbEstimate`](http://oshlacklab.com/splatter/reference/zinbEstimate.md).
For details of the ZINB-WaVE simulation see
[`zinbSimulate`](http://oshlacklab.com/splatter/reference/zinbSimulate.md).
