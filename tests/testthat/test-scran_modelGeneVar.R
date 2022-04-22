# scran_modelGeneVar.R
library(singleCellTK)
context("Testing scran_modelGeneVar.R")
data(scExample, package = "singleCellTK")

test_that(desc = "Testing scran_modelGeneVar.R", {
  sce <- runSeuratNormalizeData(sce)
  sce <- scranModelGeneVar(sce, "seuratNormData")
  
  testthat::expect_true(!is.null(rowData(sce)$scran_modelGeneVar_mean))
  testthat::expect_true(!is.null(rowData(sce)$scran_modelGeneVar_totalVariance))
  testthat::expect_true(!is.null(rowData(sce)$scran_modelGeneVar_bio))
})
