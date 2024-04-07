# Clustering Functions
library(singleCellTK)
context("Testing DEG functions")
data(scExample, package = "singleCellTK")
sce <- subsetSCECols(sce, colData = 'type != "EmptyDroplet"')
sce <- scaterlogNormCounts(sce, "logcounts")
sce <- scaterPCA(sce, useFeatureSubset = NULL)
altExp(sce, "hvg") <- sce

test_that(desc = "Testing Scran SNN with Assay", {
  sce <- runScranSNN(sce, useReducedDim = NULL, k = 8, weightType = "rank", useAssay = "logcounts",
                     clusterName = "logcounts_cluster")

  testthat::expect_true("logcounts_cluster" %in% names(colData(sce)))
})

test_that(desc = "Testing Scran SNN with PCA", {
  sce <- runScranSNN(sce, useReducedDim = "PCA",
                     clusterName = "PCA_cluster")

  testthat::expect_true("PCA_cluster" %in% names(colData(sce)))
})

test_that(desc = "Testing Scran SNN with altExp", {
  sce <- runScranSNN(sce, useReducedDim = NULL, useAltExp = "hvg", altExpAssay = "logcounts",
                     clusterName = "hvg_cluster", k = 8, weightType = "rank")

  testthat::expect_true("hvg_cluster" %in% names(colData(sce)))
})

test_that(desc = "Testing KMeans", {
  sce <- runKMeans(sce, nCenters = 2)

  testthat::expect_true("KMeans_cluster" %in% names(colData(sce)))
})
