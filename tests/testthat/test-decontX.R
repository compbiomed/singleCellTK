# decontamination algorithms
library(singleCellTK)
context("Testing decontamination algorithms")
data(scExample, package = "singleCellTK")
sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")

test_that(desc = "Testing runDecontX", {
        sceres <- runDecontX(sce)
        expect_equal(length(colData(sceres)$decontX_clusters),ncol(sce))
        expect_equal(class(colData(sceres)$decontX_contamination), "numeric")
})

test_that(desc = "runDecontX uses the decontX package, not celda", {
  sceres <- runDecontX(sce)
  runParams <- S4Vectors::metadata(sceres)$sctk$runDecontX$all_cells
  expect_equal(runParams$packageVersion,
               utils::packageDescription("decontX")$Version)
})

