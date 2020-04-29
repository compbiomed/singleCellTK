library('SummarizedExperiment')
library('SingleCellExperiment')
library('Seurat')
library('testthat')
data(sce_chcl, package = "scds")

context("Testing QC functions")

test_that("Testing runScrublet",{
  sce_testscrublet <- runScrublet(sce_chcl)
  expect_equal(class(colData(sce_testscrublet)$scrublet_score), 'numeric')
  expect_equal(class(colData(sce_testscrublet)$scrublet_call), 'logical')
  expect_equal(dim(reducedDim(sce_testscrublet,'TSNE')), c(dim(sce_testscrublet)[1],2))
  expect_equal(dim(reducedDim(sce_testscrublet,'UMAP')), c(dim(sce_testscrublet)[1],2))
})

test_that(desc = "Testing runDoubletCells", {
  sceres <- runDoubletCells(sce)
  expect_equal(length(colData(sceres)$scran_doubletCells_Score),ncol(sce))
  expect_equal(class(colData(sceres)$scran_doubletCells_Score), "numeric")
})

test_that(desc = "Testing runDoubletFinder",  {
  sceres <- runDoubletFinder(sce, seuratPcs = 1:3, seuratNfeatures = 100, seuratRes = 1)
  expect_equal(length(colData(sceres)$doubletFinder_doublet_score_Resolution_1),ncol(sce))
  expect_equal(class(colData(sceres)$doubletFinder_doublet_score_Resolution_1), "numeric")
})
