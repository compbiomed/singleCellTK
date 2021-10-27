# SingleR
library(singleCellTK)
context("Testing Cell Type Labeling functions")
data(scExample, package = "singleCellTK")
sce <- subsetSCECols(sce, colData = 'type != "EmptyDroplet"')
sce <- scaterlogNormCounts(sce, "logcounts")
rownames(sce) <- rowData(sce)$feature_name

test_that(desc = "Testing SingleR", {
  sce <- runSingleR(sce)

  testthat::expect_true("SingleR_hpca_main_pruned.labels" %in% names(colData(sce)))
  testthat::expect_true("SingleR_hpca_main_labels" %in% names(colData(sce)))
  testthat::expect_true("SingleR_hpca_main_first.labels" %in% names(colData(sce)))
  testthat::expect_true("SingleR_hpca_main_scores" %in% names(colData(sce)))
})
