# Clustering Functions
library(singleCellTK)
context("Testing DEG functions")
data(scExample, package = "singleCellTK")
sce <- subsetSCECols(sce, colData = 'type != "EmptyDroplet"')
sce <- scaterlogNormCounts(sce, "logcounts")
sce <- scaterPCA(sce, scale = TRUE)
altExp(sce, "hvg") <- sce

test_that(desc = "Testing Scran SNN with Assay", {
  sce <- runScranSNN(sce, useAssay = "logcounts",
                     clusterName = "logcounts_cluster")

  testthat::expect_true("logcounts_cluster" %in% names(colData(sce)))
})

test_that(desc = "Testing Scran SNN with PCA", {
  sce <- runScranSNN(sce, useReducedDim = "PCA",
                     clusterName = "PCA_cluster")

  testthat::expect_true("PCA_cluster" %in% names(colData(sce)))
})

test_that(desc = "Testing Scran SNN with altExp", {
  sce <- runScranSNN(sce, useAltExp = "hvg", altExpAssay = "logcounts",
                     clusterName = "hvg_cluster")

  testthat::expect_true("hvg_cluster" %in% names(colData(sce)))
})

test_that(desc = "Testing KMeans", {
  sce <- runKMeans(sce, nCenters = 2)

  testthat::expect_true("KMeans_cluster" %in% names(colData(sce)))
})
