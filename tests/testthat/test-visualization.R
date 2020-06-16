# dimensionality reduction algorithms
library(singleCellTK)
context("Testing dimensionality reduction algorithms")
data(scExample, package = "singleCellTK")
sce@assays@data$logcounts=log10(sce@assays@data$counts + 1)
sceres <- getUMAP(inSCE = sce, useAssay = "logcounts", sample = NULL, nNeighbors = 30, reducedDimName = "UMAP",
                nIterations = 20, alpha = 1, minDist = 0.01, pca = TRUE, initialDims = 50)

test_that(desc = "Testing getUMAP", {
        expect_equal(names(reducedDims(sceres)), "UMAP")
	expect_equal(nrow(reducedDim(sceres, "UMAP")), ncol(sce))
})

test_that(desc = "Testing plotSCEScatter functions", {
    p1 <- plotSCEScatter(inSCE = sceres, legendTitle = NULL,
        slot = "assays", annotation = "counts", feature = "ENSG00000251562",
        reducedDimName = "UMAP", labelClusters = FALSE)
    expect_is(p1, "ggplot")
    p2 <- plotSCEDimReduceFeatures(inSCE = sceres, feature = "ENSG00000251562",
        shape = NULL, reducedDimName = "UMAP",
        useAssay = "counts", xlab = "UMAP1", ylab = "UMAP2")
    expect_is(p2, "ggplot")
    p3 <- plotSCEDimReduceColData(inSCE = sceres, colorBy = "type",
        shape = NULL, conditionClass = "factor",
        reducedDimName = "UMAP",
        xlab = "UMAP1", ylab = "UMAP2", labelClusters = TRUE)
    expect_is(p3, "ggplot")
})

test_that(desc = "Testing plotSCEViolin functions", {
    p1 <- plotSCEViolin(inSCE = sceres, slot = "assays",
        annotation = "counts", feature = "ENSG00000251562", groupby = "type")
    expect_is(p1, "ggplot")
    p2 <- plotSCEViolinAssayData(inSCE = sceres,
        feature = "ENSG00000251562", groupby = "type")
    expect_is(p2, "ggplot")
    p3 <- plotSCEViolinColData(inSCE = sceres,
        coldata = "type", groupby = "sample")
    expect_is(p3, "ggplot")
})


