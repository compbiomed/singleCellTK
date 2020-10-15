library(singleCellTK)
context("Testing mergeSCColData")

data(scExample, package = "singleCellTK")

sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")

colData(sce)$column_name = rownames(colData(sce))
test_that(desc = "Testing mergeSCEColData", {
    sce2 <- sce
    colData(sce2)$test <- 0

    #test again
    mergedsce <- mergeSCEColData(sce, sce2)

    expect_equal(ncol(colData(sce)) + 1, ncol(colData(mergedsce)))
})



