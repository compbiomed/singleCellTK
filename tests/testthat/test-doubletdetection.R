# doublet detection algorithms
library(singleCellTK)
context("Testing doublet detection algorithms")
#sce <- mouseBrainSubsetSCE
data(sceQCExample, package = "singleCellTK")
sce <- sce[, colData(sce)$type != 'EmptyDroplet']

test_that(desc = "Testing runDoubletCells", {
	sceres <- runDoubletCells(sce)
	expect_equal(length(colData(sceres)$scran_doubletCells_Score),ncol(sce))
	expect_equal(class(colData(sceres)$scran_doubletCells_Score), "numeric")
})

test_that(desc = "Testing runDoubletFinder",  {
        sceres <- runDoubletFinder(sce, seuratPcs = 1:3, seuratNfeatures = 100, seuratRes = 1)
        expect_equal(length(colData(sceres)$doubletFinder_doublet_score_Resolution_1),ncol(sce))
        expect_equal(class(colData(sceres)$doubletFinder_doublet_score_Resolution_1), "numeric")
})

