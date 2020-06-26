library('SummarizedExperiment')
library('SingleCellExperiment')
library('singleCellTK')
library('Seurat')
library('testthat')

data(scExample, package = "singleCellTK")
fullsce <- sce  ## sce including empty droplet
sce <- sce[, colData(sce)$type != 'EmptyDroplet']

context("Testing QC functions")

test_that("Testing scds",{
  sce <- runCxdsBcdsHybrid(sce, estNdbl = TRUE)
  expect_equal(class(colData(sce)$scds_hybrid_score), 'numeric')
  #expect_equal(class(colData(sce)$doublet_true_labels), 'character')
  expect_equal(class(colData(sce)$scds_hybrid_call), 'logical')
})

test_that("Testing emptydrops",{
  fullsce <- runEmptyDrops(fullsce)
  #expect_equal(class(colData(sce)$doublet_true_labels), 'character')
  expect_equal(class(colData(fullsce)$dropletUtils_emptyDrops_total), 'integer')
  expect_equal(class(colData(fullsce)$dropletUtils_emptyDrops_logprob), 'numeric')
  expect_equal(class(colData(fullsce)$dropletUtils_emptyDrops_pvalue), 'numeric')
  expect_equal(class(colData(fullsce)$dropletUtils_emptyDrops_limited), 'logical')
  expect_equal(class(colData(fullsce)$dropletUtils_emptyDrops_fdr), 'numeric')
})


test_that(desc = "Testing scran", {
  sce <- runDoubletCells(sce)
  expect_equal(length(colData(sce)$scran_doubletCells_Score),ncol(sce))
  expect_equal(class(colData(sce)$scran_doubletCells_Score), "numeric")
}) 

test_that(desc = "Testing DoubletFinder",  {
  sce <- runDoubletFinder(sce, seuratPcs = 1:3, seuratNfeatures = 300, seuratRes = 1,
	 verbose = FALSE, seed = 12345)
  expect_equal(length(colData(sce)$doubletFinder_doublet_score_Resolution_1),ncol(sce))
  expect_equal(class(colData(sce)$doubletFinder_doublet_score_Resolution_1), "numeric")
})


test_that(desc = "Testing runDoubletCells", {
  sceres <- runDoubletCells(sce)
  expect_equal(length(colData(sceres)$scran_doubletCells_Score),ncol(sce))
  expect_equal(class(colData(sceres)$scran_doubletCells_Score), "numeric")
})

test_that("Testing scrublet",{
  if (!reticulate::py_module_available("scanpy") || (!reticulate::py_module_available("scrublet"))){
    skip("scrublet or scanpy not available. Skipping testing importOptimus")
  }
  sce <- runScrublet(sce)
  expect_equal(class(colData(sce)$scrublet_score), 'numeric')
  expect_equal(class(colData(sce)$scrublet_call), 'logical')
  expect_equal(dim(reducedDim(sce,'scrublet_TSNE')), c(ncol(sce),2))
  expect_equal(dim(reducedDim(sce,'scrublet_UMAP')), c(ncol(sce),2))
})
