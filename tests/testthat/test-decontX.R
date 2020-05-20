# decontamination algorithms
library(singleCellTK)
context("Testing decontamination algorithms")
#sce <- mouseBrainSubsetSCE
data(sceQCExample, package = "singleCellTK")
sce <- sce[, colData(sce)$type != 'EmptyDroplet']

test_that(desc = "Testing runDecontX", {
        sceres <- runDecontX(sce)
        expect_equal(length(colData(sceres)$decontX_clusters),ncol(sce))
        expect_equal(class(colData(sceres)$decontX_contamination), "numeric")
})

