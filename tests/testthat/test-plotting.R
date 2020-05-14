library(singleCellTK)
context("Testing plotting functions")

test_that(desc = "Testing plotSCEScatter functions", {
    p1 <- plotSCEScatter(inSCE = mouseBrainSubsetSCE, legendTitle = NULL,
        slot = "assays", annotation = "counts", feature = "Tspan12",
        reducedDimName = "TSNE_counts", labelClusters = FALSE)
    expect_is(p1, "ggplot")
    p2 <- plotSCEDimReduceFeatures(inSCE = mouseBrainSubsetSCE, feature = "Sox2",
        shape = "No Shape", reducedDimName = "TSNE_counts",
        useAssay = "counts", xlab = "tSNE1", ylab = "tSNE2")
    expect_is(p2, "ggplot")
    p3 <- plotSCEDimReduceColData(inSCE = mouseBrainSubsetSCE, colorBy = "tissue",
        shape = "No Shape", conditionClass = "factor",
        reducedDimName = "TSNE_counts",
        xlab = "tSNE1", ylab = "tSNE2", labelClusters = TRUE)
    expect_is(p3, "ggplot")
})

test_that(desc = "Testing plotSCEScatter functions", {
    p1 <- plotSCEViolin(inSCE = mouseBrainSubsetSCE, slot = "assays",
        annotation = "counts", feature = "Sox2", groupby = "sex")
    expect_is(p1, "ggplot")
    p2 <- plotSCEViolinAssayData(inSCE = mouseBrainSubsetSCE,
        feature = "Sox2", groupby = "sex")
    expect_is(p2, "ggplot")
    p3 <- plotSCEViolinColData(inSCE = mouseBrainSubsetSCE,
        coldata = "age", groupby = "sex")
    expect_is(p3, "ggplot")
})

