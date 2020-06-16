
# QC algorithms for all barcodes, pre-filter
library(singleCellTK)
context("Testing dropletUtils algorithms")
data(scExample, package = "singleCellTK")

test_that(desc = "Testing runBarcodeRankDrops", {
        sceres <- runBarcodeRankDrops(inSCE = sce)
        expect_equal(length(colData(sceres)$dropletUtils_BarcodeRank_Inflection),ncol(sce))
        expect_equal(class(colData(sceres)$dropletUtils_BarcodeRank_Inflection), "integer")
})

