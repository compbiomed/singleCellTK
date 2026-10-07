# Trajectory Analysis
library(singleCellTK)
context("Testing Trajectory analysis function")

data("scExample", package = "singleCellTK")
sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
sce <- runNormalization(sce,
                        normalizationMethod = "LogNormalize",
                        useAssay = "counts",
                        outAssayName = "logcounts")
sce <- runDimReduce(inSCE = sce,
                    method = "scaterPCA",
                    useAssay = "logcounts",
                    reducedDimName = "PCA")
sce <- runDimReduce(inSCE = sce,
                    method = "rTSNE",
                    useReducedDim = "PCA",
                    reducedDimName = "TSNE")

test_that(desc = "Testing TSCAN", {
  sce <- runTSCAN(inSCE = sce, useReducedDim = "PCA", seed = NULL)
  terminalNodes <- listTSCANTerminalNodes(sce)
  sce <- runTSCANDEG(inSCE = sce, pathIndex = terminalNodes[1])
  sce <- runTSCANClusterDEAnalysis(inSCE = sce, useCluster = 1)
  testthat::expect_true(!is.null(names(
                          getTSCANResults(sce,
                                          analysisName = "Pseudotime"))))
  testthat::expect_true(!is.null(names(
                          getTSCANResults(sce,
                                          analysisName = "DEG",
                                          terminalNodes[1]))))
  testthat::expect_true(!is.null(names(
                          getTSCANResults(sce,
                                          analysisName = "ClusterDEAnalysis",
                                          1))))


  TSCANResultsPlot <- plotTSCANResults(inSCE = sce, useReducedDim = "TSNE")
  TSCANPseudotimeHeatmapPlot1 <- plotTSCANPseudotimeHeatmap(inSCE = sce,
                                                           pathIndex = terminalNodes[1],
                                                           direction = "both",
                                                           log2fcThreshold = 0.001,
                                                           topN = 5)
  TSCANPseudotimeHeatmapPlot2 <- plotTSCANPseudotimeHeatmap(inSCE = sce,
                                                            pathIndex = terminalNodes[1],
                                                            direction = "inc",
                                                            log2fcThreshold = 0.001,
                                                            topN = 5)
  TSCANPseudotimeGenesPlot <- plotTSCANPseudotimeGenes(inSCE = sce,
                                                       pathIndex = terminalNodes[1],
                                                       direction = "increasing")
  ClusterPseudoPlot <- plotTSCANClusterPseudo(inSCE = sce,
                                              useCluster = 1,
                                              useReducedDim = "TSNE")
  TSCANBranchGenesPlot <- plotTSCANClusterDEG(inSCE = sce,
                                             useCluster = 1,
                                             useReducedDim = "TSNE")
  TSCANFeaturePlot <- plotTSCANDimReduceFeatures(inSCE = sce,
                                                 feature = "CD74",
                                                 by = "feature_name",
                                                 useReducedDim = "TSNE",
                                                 featureDisplay = "feature_name")

  testthat::expect_is(TSCANResultsPlot, "ggplot")
  testthat::expect_is(TSCANPseudotimeHeatmapPlot1, "Heatmap")
  testthat::expect_is(TSCANPseudotimeHeatmapPlot2, "Heatmap")
  testthat::expect_is(TSCANPseudotimeGenesPlot, "ggplot")
  testthat::expect_is(ClusterPseudoPlot, "ggplot")
  testthat::expect_is(TSCANFeaturePlot, "ggplot")
})



test_that(desc = "TSCAN cluster centroids match base R without deprecations", {
  sceT <- suppressWarnings(runTSCAN(inSCE = sce, useReducedDim = "PCA",
                                    seed = NULL))
  expect_no_warning(plotTSCANResults(inSCE = sceT, useReducedDim = "TSNE"),
                    class = "deprecatedWarning")
  clusters <- sceT$TSCAN_clusters
  emb <- reducedDim(sceT, "PCA")
  ids <- levels(factor(clusters))
  expected <- t(vapply(ids,
                       function(cl) colMeans(emb[clusters == cl, , drop = FALSE]),
                       numeric(ncol(emb))))
  expect_equal(.clusterCentroids(sceT, clusters, "PCA"), expected)
})

test_that(desc = "TSCAN cluster centroids ignore cells without a cluster", {
  clusters <- rep(c("a", "b", NA), length.out = ncol(sce))
  emb <- reducedDim(sce, "PCA")
  expected <- rbind(a = colMeans(emb[which(clusters == "a"), ]),
                    b = colMeans(emb[which(clusters == "b"), ]))
  expect_equal(.clusterCentroids(sce, clusters, "PCA"), expected)
})
