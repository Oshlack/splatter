library(splatter)
library(scran)
library(scater)

params = newSplatParams()
params = setParams(params,
                   nGenes = 500,
                   batchCells = rep(1000, 2),
                   batch.facLoc = 0.05,
                   batch.facScale = 0.05,
                   batch.rmEffect = FALSE,
                   splits.per.level = c(4, 2, 2),
                   de.prob.per.level = c(0.05, 0.0333, 0.0167),
                   de.facLoc.per.level = 0.1,
                   de.facScale.per.level = 0.4,
                   seed = 1)

data = splatSimulate(params, method = "hierarchical", verbose = FALSE)

data = logNormCounts(data)
data = runPCA(data)
data = runUMAP(data)

plot(hclust(dist(t(as.matrix(rowData(data)[, grepl("DEFacGroup", colnames(rowData(data))), drop = FALSE])))))

means = sapply(sort(unique(data$Group)), function(x) rowMeans(logcounts(data)[, data$Group == x]))
colnames(means) = sort(unique(data$Group))
plot(hclust(dist(t(means))))

cowplot::plot_grid(plotUMAP(data, colour_by = "Group"), plotUMAP(data, colour_by = "Batch"), ncol = 2)
