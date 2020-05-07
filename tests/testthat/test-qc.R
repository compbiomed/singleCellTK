library('SummarizedExperiment')
library('SingleCellExperiment')
library('singleCellTK')
library('Seurat')
library('testthat')

data(sce_chcl, package = "scds")

context("Testing QC functions")

test_that("Testing scrublet",{
  sce <- runScrublet(sce_chcl)
  expect_equal(class(colData(sce)$scrublet_score), 'numeric')
  expect_equal(class(colData(sce)$scrublet_call), 'logical')
  expect_equal(dim(reducedDim(sce,'TSNE')), c(dim(sce)[1],2))
  expect_equal(dim(reducedDim(sce,'UMAP')), c(dim(sce)[1],2))
})

test_that("Testing scds",{
  sce <- runCxdsBcdsHybrid(sce_chcl)
  expect_equal(class(colData(sce)$scds_hybrid_score), 'numeric')
  expect_equal(class(colData(sce)$doublet_true_labels), 'character')
})

test_that("Testing emptydrops",{
  sce <- runEmptyDrops(sce_chcl)
  expect_equal(class(colData(sce)$doublet_true_labels), 'character')
  expect_equal(class(colData(sce)$dropletUtils_emptyDrops_total), 'integer')
  expect_equal(class(colData(sce)$dropletUtils_emptyDrops_logprob), 'numeric')
  expect_equal(class(colData(sce)$dropletUtils_emptyDrops_pvalue), 'numeric')
  expect_equal(class(colData(sce)$dropletUtils_emptyDrops_limited), 'logical')
  expect_equal(class(colData(sce)$dropletUtils_emptyDrops_fdr), 'numeric')
})


test_that(desc = "Testing scran", {
  sce <- runDoubletCells(sce_chcl)
  expect_equal(length(colData(sce)$scran_doubletCells_Score),ncol(sce_chcl))
  expect_equal(class(colData(sce)$scran_doubletCells_Score), "numeric")
}) 

test_that(desc = "Testing DoubletFinder",  {
  sce <- runDoubletFinder(sce_chcl, seuratPcs = 1:3, seuratNfeatures = 100, seuratRes = 1)
  expect_equal(length(colData(sce)$doubletFinder_doublet_score_Resolution_1),ncol(sce_chcl))
  expect_equal(class(colData(sce)$doubletFinder_doublet_score_Resolution_1), "numeric")
})


