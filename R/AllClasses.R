setClass("Params",
         contains = "VIRTUAL",
         slots = c(nGenes = "numeric",
                   nCells = "numeric",
                   seed = "numeric"),
         prototype = prototype(nGenes = 10000, nCells = 100,
                               seed = sample(1:1e6, 1)))

setClass("SimpleParams",
         contains = "Params",
         slots = c(mean.shape = "numeric",
                   mean.rate = "numeric",
                   count.disp = "numeric"),
         prototype = prototype(mean.shape = 0.4, mean.rate = 0.3,
                               count.disp = 0.1))