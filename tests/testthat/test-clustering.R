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

samePartition <- function(a, b) {
  length(unique(paste(a, b))) == length(unique(a)) &&
    length(unique(a)) == length(unique(b))
}

test_that(desc = "Scran SNN builds graphs without deprecations", {
  expect_no_warning(
    res <- runScranSNN(sce, useReducedDim = "PCA", k = 8, nComp = 5,
                       weightType = "jaccard", algorithm = "walktrap",
                       clusterName = "pca_snn"),
    class = "deprecatedWarning"
  )
  pcs <- reducedDim(sce, "PCA")[, seq(5)]
  g <- bluster::makeSNNGraph(pcs, k = 8, type = "jaccard")
  expected <- igraph::cluster_walktrap(g)$membership
  expect_true(samePartition(as.integer(res$pca_snn), expected))
  expect_no_warning(
    runScranSNN(sce, useReducedDim = NULL, useAssay = "logcounts", k = 8,
                clusterName = "assay_snn"),
    class = "deprecatedWarning"
  )
  expect_no_warning(
    runScranSNN(sce, useReducedDim = NULL, useAltExp = "hvg",
                altExpAssay = "logcounts", k = 8, clusterName = "ae_snn"),
    class = "deprecatedWarning"
  )
})

test_that(desc = "Scran SNN uses the reducedDim of an altExp", {
  res <- suppressWarnings(
    runScranSNN(sce, useReducedDim = NULL, useAltExp = "hvg",
                altExpRedDim = "PCA", k = 8, clusterName = "ae_pca_snn"))
  expect_true("ae_pca_snn" %in% names(colData(res)))
})
