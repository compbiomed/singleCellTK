# Trajectory Analysis
library(singleCellTK)
context("Testing Trajectory analysis function")

data("scExample", package = "singleCellTK")
sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
rowData(sce)$Symbol <- rowData(sce)$feature_name
rownames(sce) <- rowData(sce)$Symbol
sce <- scaterlogNormCounts(sce, assayName = "logcounts")
sce <- runDimReduce(inSCE = sce, method = "scaterPCA", useAssay = "logcounts", reducedDimName = "PCA")
sce <- runDimReduce(inSCE = sce, method = "rTSNE", useReducedDim = "PCA", reducedDimName = "TSNE")

test_that(desc = "Testing TSCAN", {
  sce <- runTSCAN (inSCE = sce, useReducedDim = "PCA", seed = NULL)
  sce <- runTSCANDEG(inSCE = sce, pathIndex = 4)
  sce <- runTSCANClusterDEAnalysis(inSCE = sce, useClusters = 5)
  
  testthat::expect_true(!is.null(names(getTSCANResults(sce, analysisName = "Pseudotime"))))
  testthat::expect_true(!is.null(names(getTSCANResults(sce, analysisName = "DEG"))))
  testthat::expect_true(!is.null(names(getTSCANResults(sce, analysisName = "ClusterDEAnalysis"))))
  
  
  TSCANResultsPlot <- plotTSCANResults(inSCE = sce, useReducedDim = "TSNE")
  TSCANPseudotimeHeatmapPlot <- plotTSCANPseudotimeHeatmap(inSCE = sce, pathIndex = 4,topN = 5)
  TSCANPseudotimeGenesPlot <- plotTSCANPseudotimeGenes(inSCE = sce, pathIndex = 4, direction = "increasing")
  ClusterPseudoPlot <- plotClusterPseudo(inSCE = sce, useClusters = 5, pathIndex = NULL, useReducedDim = "TSNE")
  TSCANDEgenesPlot <- plotTSCANDEgenes(inSCE = sce, geneSymbol = "CD74", useReducedDim = "TSNE")
  
  testthat::expect_is(TSCANResultsPlot, "ggplot")
  testthat::expect_is(TSCANPseudotimeHeatmapPlot, "ggplot")
  testthat::expect_is(TSCANPseudotimeGenesPlot, "ggplot")
  testthat::expect_is(ClusterPseudoPlot, "ggplot")
  testthat::expect_is(TSCANDEgenesPlot, "ggplot")
  
  
})


