# seurat functions
library(singleCellTK)
context("Testing seurat functions")
data(scExample, package = "singleCellTK")

test_that(desc = "Testing standard seurat workflow", {
  sce <- seuratNormalizeData(sce)
  testthat::expect_true("seuratNormData" %in% assayNames(sce))
  
  sce <- seuratScaleData(sce)
  testthat::expect_true("seuratScaledData" %in% assayNames(sce))
  
  sce <- seuratFindHVG(sce)
  testthat::expect_true("seurat_variableFeatures_vst_varianceStandardized"
                        %in% names(rowData(sce)))
})


