library('SummarizedExperiment')
library('SingleCellExperiment')
library('Seurat')
library('testthat')
data(sce_chcl, package = "scds")

context("Testing QC functions")

test_that("Testing runScrublet",{
  sce <- runScrublet(sce_chcl)
  expect_equal(class(colData(sce)$scrublet_score), 'numeric')
  expect_equal(class(colData(sce)$scrublet_call), 'logical')
  expect_equal(dim(reducedDim(sce,'TSNE')), c(dim(sce)[1],2))
  expect_equal(dim(reducedDim(sce,'UMAP')), c(dim(sce)[1],2))
})

test_that("Testing runCxdsBcdsHybrid",{
  sce <- runCxdsBcdsHybrid(sce_chcl)
  expect_equal(class(colData(sce)$scds_hybrid_score), 'numeric')
  expect_equal(class(colData(sce)$doublet_true_labels), 'character')
})

test_that(desc = "Testing runDoubletCells", {
  sce <- runDoubletCells(sce_chcl)
  expect_equal(length(colData(sce)$scran_doubletCells_Score),ncol(sce_chcl))
  expect_equal(class(colData(sce)$scran_doubletCells_Score), "numeric")
})

test_that(desc = "Testing runDoubletFinder",  {
  sce <- runDoubletFinder(sce_chcl, seuratPcs = 1:3, seuratNfeatures = 100, seuratRes = 1)
  expect_equal(length(colData(sce)$doubletFinder_doublet_score_Resolution_1),ncol(sce_chcl))
  expect_equal(class(colData(sce)$doubletFinder_doublet_score_Resolution_1), "numeric")
})


