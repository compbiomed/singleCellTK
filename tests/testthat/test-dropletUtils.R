
# QC algorithms for all barcodes, pre-filter
library(singleCellTK)
context("Testing dropletUtils algorithms")
data(scExample, package = "singleCellTK")

test_that(desc = "Testing runBarcodeRankDrops", {
        sce <- runBarcodeRankDrops(inSCE = sce)
        expect_equal(length(colData(sce)$dropletUtils_BarcodeRank_Inflection),ncol(sce))
        expect_equal(class(colData(sce)$dropletUtils_BarcodeRank_Inflection), "integer")
        gg <- plotBarcodeRankDropsResults(sce)
        expect_is(gg$scatterBarcodeRank, "ggplot")
})

