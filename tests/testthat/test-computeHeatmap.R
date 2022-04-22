# computeHeatmap.R
library(singleCellTK)
library(ggplot2)
context("Testing computeHeatmap.R")
data(scExample, package = "singleCellTK")

test_that(desc = "Testing computeHeatmap.R", {
  sce <- runSeuratNormalizeData(sce)
  sce <- runSeuratScaleData(sce)
  sce <- runSeuratFindHVG(sce)
  sce <- runSeuratPCA(sce)
  heatmap_gTree_objects <- computeHeatmap(sce, useAssay = "seuratScaledData", dims = 10)
  heatmapPlot <- singleCellTK:::.plotHeatmapMulti(heatmap_gTree_objects)
  
  testthat::expect_equal(length(heatmap_gTree_objects), 10)
  testthat::expect_true(is.ggplot(heatmapPlot))
})
