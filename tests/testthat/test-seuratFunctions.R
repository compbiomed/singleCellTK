# seurat functions
library(singleCellTK)
library(ggplot2)
context("Testing seurat functions")
data(scExample, package = "singleCellTK")

test_that(desc = "Testing standard seurat workflow", {
  # Test Normalization Method
  sce <- runSeuratNormalizeData(sce)
  testthat::expect_true("seuratNormData" %in% assayNames(sce))
  
  # Test Scaling Method
  sce <- runSeuratScaleData(sce)
  testthat::expect_true("seuratScaledData" %in% assayNames(sce))
  
  # Test Feature Selection Method
  sce <- runSeuratFindHVG(sce)
  testthat::expect_true("seurat_variableFeatures_vst_varianceStandardized"
                        %in% names(rowData(sce)))
  
  # Test PCA
  sce <- runSeuratPCA(sce, useAssay = "seuratScaledData")
  testthat::expect_true("seuratPCA" %in% reducedDimNames(sce))
  
  # Test ICA
  sce <- runSeuratICA(sce, useAssay = "seuratScaledData")
  testthat::expect_true("seuratICA" %in% reducedDimNames(sce))
  
  # Test JackStraw on PCA
  sce <- suppressWarnings({runSeuratJackStraw(sce, useAssay = "seuratScaledData")})
  testthat::expect_true(!is.null(metadata(sce)$seurat$obj@reductions$pca@jackstraw))
  
  jsPlot <- plotSeuratJackStraw(sce)
  testthat::expect_true(is.ggplot(jsPlot))
  
  # Test HVG Plot 
  hvgPlot <- plotSeuratHVG(sce, labelPoints = 10)
  testthat::expect_true(is.ggplot(hvgPlot))
  
  # Test PCA Component Plot
  pcaPlot <- plotSeuratReduction(sce, useReduction = "pca")
  testthat::expect_true(is.ggplot(pcaPlot))
  pcaPlot <- plotSeuratReduction(sce, useReduction = "pca", showLegend = TRUE)
  testthat::expect_true(is.ggplot(pcaPlot))
  
  # Test Clustering
  sce <- runSeuratFindClusters(sce, useAssay = "seuratScaledData")
  testthat::expect_true("Seurat_louvain_Resolution0.8" %in% names(colData(sce)))
  
  # Test TSNE
  sce <- runSeuratTSNE(sce, useReduction = "pca")
  testthat::expect_true("seuratTSNE" %in% reducedDimNames(sce))
  
  # Test UMAP
  sce <- suppressWarnings({runSeuratUMAP(sce, useReduction = "pca")})
  testthat::expect_true("seuratUMAP" %in% reducedDimNames(sce))
  
  # Test ElbowPlot on PCA
  elbowPlot <- plotSeuratElbow(sce)
  testthat::expect_true(inherits(elbowPlot, "plotly"))
  
  # Test Heatmap on PCA
  heatmapPlot <- runSeuratHeatmap(sce, useAssay = "seuratScaledData", dims = 4, fast = FALSE)
  testthat::expect_true(is.ggplot(heatmapPlot))
  
  # Test SCTransform (requires removal of zero rowSums/colSums)
  zeroCols <- which(colSums(assay(sce, "counts")) == 0)
  sce2 <- sce[, -zeroCols]
  metadata(sce2)$seurat <- NULL
  sce2 <- runSeuratSCTransform(sce2)
  testthat::expect_true("SCTCounts" %in% assayNames(sce2))
  
  # Test Batch-Correction
  sce <- suppressWarnings(runSeuratIntegration(sce, batch = "type", kAnchor = 10, kFilter = 4, kWeight = 5, ndims = 10))
  testthat::expect_true("SeuratIntegratedAssay" %in% altExpNames(sce))
  
  # Test Marker Selection - Conserved Between 2 Groups
  cells1 <- colnames(sce)[which(colData(sce)$type == "Singlet")]
  cells2 <- colnames(sce)[which(colData(sce)$type == "Doublet")]
  sce <- runSeuratFindMarkers(sce, cells1 = cells1, cells2 = cells2, group1 = "Singlet", group2 = "Doublet", conserved = TRUE)
  testthat::expect_true("seuratMarkers" %in% names(metadata(sce)))
  
  # Test Marker Selection - Standard Between All Groups
  metadata(sce)$seurat
  sce <- runSeuratFindMarkers(sce, allGroup = "type")
  testthat::expect_true("seuratMarkers" %in% names(metadata(sce)))
  
  # Test Plots for Top Marker Genes
  top2Features <- metadata(sce)$seuratMarkers[1:2, ]$gene.id
  genePlot1 <- plotSeuratGenes(sce, plotType = "ridge", features = top2Features, groupVariable = "type", combine = TRUE)
  genePlot2 <- plotSeuratGenes(sce, plotType = "violin", features = top2Features, groupVariable = "type", combine = TRUE)
  genePlot3 <- plotSeuratGenes(sce, plotType = "heatmap", features = top2Features, groupVariable = "type")
  genePlot4 <- plotSeuratGenes(sce, plotType = "feature", features = top2Features, groupVariable = "type", combine = TRUE)
  genePlot5 <- plotSeuratGenes(sce, plotType = "dot", features = top2Features, groupVariable = "type")
  testthat::expect_true(is.ggplot(genePlot1))
  testthat::expect_true(is.ggplot(genePlot2))
  testthat::expect_true(is.ggplot(genePlot3))
  testthat::expect_true(is.ggplot(genePlot4))
  testthat::expect_true(is.ggplot(genePlot5))
})


