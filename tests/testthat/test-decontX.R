# decontamination algorithms
library(singleCellTK)
context("Testing decontamination algorithms")
sce <- emptyDropsSceExample

test_that(desc = "Testing runDoubletCells", {
        sceres <- runDecontX(sce)
        expect_equal(length(colData(sceres)$decontX_clusters),ncol(sce))
        expect_equal(class(colData(sceres)$decontX_contamination), "numeric")
})

