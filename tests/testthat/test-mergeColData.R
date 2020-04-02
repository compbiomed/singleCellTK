library(singleCellTK)
context("Testing fxn")
sce <- mouseBrainSubsetSCE

test_that(desc = "Testing runDoubletCells", {
    sce2 <- sce
    colData(sce2)$test <- 0

    #test again
    mergedsce <- mergeSCEColData(sce, sce2)

    expect_equal(ncol(colData(sce) + 1), ncol(colData(mergedsce)))
}

