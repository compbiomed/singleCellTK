# seurat functions
library(singleCellTK)
library(ggplot2)
context("Testing seurat functions")
data(scExample, package = "singleCellTK")

test_that(desc = "Testing standard seurat workflow", {
  # Test Normalization Method
  sce <- seuratNormalizeData(sce)
  testthat::expect_true("seuratNormData" %in% assayNames(sce))
  
  # Test Scaling Method
  sce <- seuratScaleData(sce)
  testthat::expect_true("seuratScaledData" %in% assayNames(sce))
  
  # Test Feature Selection Method
  sce <- seuratFindHVG(sce)
  testthat::expect_true("seurat_variableFeatures_vst_varianceStandardized"
                        %in% names(rowData(sce)))
  
  # Test PCA
  sce <- seuratPCA(sce, useAssay = "seuratScaledData")
  testthat::expect_true("seuratPCA" %in% reducedDimNames(sce))
  
  # Test ICA
  sce <- seuratICA(sce, useAssay = "seuratScaledData")
  testthat::expect_true("seuratICA" %in% reducedDimNames(sce))
  
  # Test JackStraw on PCA
  sce <- suppressWarnings({seuratComputeJackStraw(sce, useAssay = "seuratScaledData")})
  testthat::expect_true(!is.null(metadata(sce)$seurat$obj@reductions$pca@jackstraw))
  
  jsPlot <- seuratJackStrawPlot(sce)
  testthat::expect_true(is.ggplot(jsPlot))
  
  # Test HVG Plot 
  hvgPlot <- seuratPlotHVG(sce, labelPoints = 10)
  testthat::expect_true(is.ggplot(hvgPlot))
  
  # Test PCA Component Plot
  pcaPlot <- seuratReductionPlot(sce, useReduction = "pca")
  testthat::expect_true(is.ggplot(pcaPlot))
  pcaPlot <- seuratReductionPlot(sce, useReduction = "pca", showLegend = TRUE)
  testthat::expect_true(is.ggplot(pcaPlot))
  
  # Test Clustering
  sce <- seuratFindClusters(sce, useAssay = "seuratScaledData")
  testthat::expect_true("Seurat_louvain_Resolution0.8" %in% names(colData(sce)))
  
  # Test TSNE
  sce <- seuratRunTSNE(sce, useReduction = "pca")
  testthat::expect_true("seuratTSNE" %in% reducedDimNames(sce))
  
  # Test UMAP
  sce <- suppressWarnings({seuratRunUMAP(sce, useReduction = "pca")})
  testthat::expect_true("seuratUMAP" %in% reducedDimNames(sce))
  
  # Test ElbowPlot on PCA
  elbowPlot <- seuratElbowPlot(sce)
  testthat::expect_true(inherits(elbowPlot, "plotly"))
  
  # Test Heatmap on PCA
  heatmapPlot <- seuratComputeHeatmap(sce, useAssay = "seuratScaledData", dims = 4, fast = FALSE)
  testthat::expect_true(is.ggplot(heatmapPlot))
  
  # Test SCTransform (requires removal of zero rowSums/colSums)
  zeroCols <- which(colSums(assay(sce, "counts")) == 0)
  sce2 <- sce[, -zeroCols]
  metadata(sce2)$seurat <- NULL
  sce2 <- seuratSCTransform(sce2)
  testthat::expect_true("SCTCounts" %in% assayNames(sce2))
  
  # Test Batch-Correction
  sce <- suppressWarnings(seuratIntegration(sce, batch = "type", kAnchor = 10, kFilter = 4, kWeight = 5, ndims = 10))
  testthat::expect_true("SeuratIntegratedAssay" %in% altExpNames(sce))
  
  # Test Marker Selection - Conserved Between 2 Groups
  cells1 <- colnames(sce)[which(colData(sce)$type == "Singlet")]
  cells2 <- colnames(sce)[which(colData(sce)$type == "Doublet")]
  sce <- seuratFindMarkers(sce, cells1 = cells1, cells2 = cells2, group1 = "Singlet", group2 = "Doublet", conserved = TRUE)
  testthat::expect_true("seuratMarkers" %in% names(metadata(sce)))
  
  # Test Marker Selection - Standard Between All Groups
  metadata(sce)$seurat
  sce <- seuratFindMarkers(sce, allGroup = "type")
  testthat::expect_true("seuratMarkers" %in% names(metadata(sce)))
  
  # Test Plots for Top Marker Genes
  top2Features <- metadata(sce)$seuratMarkers[1:2, ]$gene.id
  genePlot1 <- seuratGenePlot(sce, plotType = "ridge", features = top2Features, groupVariable = "type")
  genePlot2 <- seuratGenePlot(sce, plotType = "violin", features = top2Features, groupVariable = "type")
  genePlot3 <- seuratGenePlot(sce, plotType = "heatmap", features = top2Features, groupVariable = "type")
  genePlot4 <- seuratGenePlot(sce, plotType = "feature", features = top2Features, groupVariable = "type")
  genePlot5 <- seuratGenePlot(sce, plotType = "dot", features = top2Features, groupVariable = "type")
  testthat::expect_true(is.ggplot(genePlot1))
  testthat::expect_true(is.ggplot(genePlot2))
  testthat::expect_true(is.ggplot(genePlot3))
  testthat::expect_true(is.ggplot(genePlot4))
  testthat::expect_true(is.ggplot(genePlot5))
})


